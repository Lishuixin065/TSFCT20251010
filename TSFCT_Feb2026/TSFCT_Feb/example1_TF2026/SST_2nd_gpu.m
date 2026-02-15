function [Ts, tfr, tfrtic, omega] = SST_2nd_gpu(x, Hz, s)
    %% GPU版本 - 完全向量化，所有步骤在GPU上
    [xrow, ~] = size(x);
    if (xrow ~= 1)
        x = x';
    end
    N = length(x);

    % GPU检查
    if gpuDeviceCount == 0
        error('GPU not available - this version requires GPU');
    end

    % 选择GPU设备
    gpuDevice(1);

    %% 基本参数
    if mod(N, 2) == 0
        f = (0:N/2) * Hz / N;
    else
        f = (0:(N-1)/2) * Hz / N;
    end
    tfrtic = f;

    Lh = N;
    ht = -N:N;
    time = ht / Hz;
    tLen = N;
    fLen = length(tfrtic);

    %% 数据转移到GPU
    x_gpu = gpuArray(single(x));
    time_gpu = gpuArray(single(time));

    %% 预计算窗口函数
    s_sq = single(s^2);
    factor = single(-2 * pi^2 * s_sq);
    
    time_sq = time_gpu.^2;
    exp_term = exp(factor * time_sq);
    g_all = exp_term;
    tg_all = exp_term .* time_gpu;
    ttg_all = exp_term .* time_sq;

    %% STFT计算
    tfr = gpuArray.zeros(fLen, tLen, 'single');
    tfr0 = gpuArray.zeros(fLen, tLen, 'single');
    tfr1 = gpuArray.zeros(fLen, tLen, 'single');
    tfr2 = gpuArray.zeros(fLen, tLen, 'single');
    
    % 预分配临时数组
    tf0 = gpuArray.zeros(N, 1, 'single');
    tf1 = gpuArray.zeros(N, 1, 'single');
    tf2 = gpuArray.zeros(N, 1, 'single');
    
    for tidx = 1:tLen
        ti = tidx;
        
        % 计算tau范围
        tau = -min([round(N/2)-1, ti-1]):min([round(N/2)-1, N-ti]);
        
        if ~isempty(tau)
            % 计算索引
            indices = rem(N + tau, N) + 1;
            
            % 提取信号和窗口函数值
            x_seg = x_gpu(ti + tau);
            g_vals = g_all(Lh + 1 + tau);
            tg_vals = tg_all(Lh + 1 + tau);
            ttg_vals = ttg_all(Lh + 1 + tau);
            
            % 重置临时数组
            tf0(:) = 0;
            tf1(:) = 0;
            tf2(:) = 0;
            
            % 赋值
            tf0(indices) = x_seg .* g_vals;
            tf1(indices) = x_seg .* tg_vals;
            tf2(indices) = x_seg .* ttg_vals;
            
            % FFT
            TFR0 = fft(tf0);
            TFR1 = fft(tf1);
            TFR2 = fft(tf2);
            
            % 存储结果
            tfr(:, tidx) = TFR0(1:fLen);
            tfr0(:, tidx) = TFR0(1:fLen);
            tfr1(:, tidx) = TFR1(1:fLen);
            tfr2(:, tidx) = TFR2(1:fLen);
        end
    end

    %% 频率重分配
    del_f = single(tfrtic(2) - tfrtic(1));
    del_f_gpu = gpuArray(del_f);
    
    % 计算分母，避免除零
    denominator = tfr1 .* tfr1 - tfr0 .* tfr2;
    eps_val = eps('single');
    denominator = denominator + (denominator == 0) * eps_val;
    
    % 计算瞬时频率
    temp = (tfr0 .* tfr1) ./ denominator;
    omega1 = real(temp ./ (single(2.0 * pi * 1i)));
    omega1 = round(omega1 ./ del_f_gpu);
    
    % 构建omega矩阵
    eta_vec = (1:fLen)';
    eta_gpu = gpuArray(single(eta_vec));
    omega = bsxfun(@minus, eta_gpu, omega1);

    %% 同步压缩 - 完全向量化GPU实现
    Ts = gpuArray.zeros(fLen, tLen, 'single');
    
    % 获取所有有效点的线性索引
    valid_mask = abs(tfr) > eps_val;
    
    if any(valid_mask(:))
        % 获取所有有效点的坐标
        [eta_indices, time_indices] = find(valid_mask);
        
        % 获取对应的tfr值和omega值
        linear_idx = eta_indices + (time_indices - 1) * fLen;
        tfr_vals = tfr(linear_idx);
        omega_vals = omega(linear_idx);
        
        % 计算目标行
        target_rows = round(omega_vals);
        
        % 筛选有效目标
        valid_targets = (target_rows >= 1) & (target_rows <= fLen);
        
        if any(valid_targets)
           
            time_valid = time_indices(valid_targets);
            target_rows = target_rows(valid_targets);
            tfr_vals = tfr_vals(valid_targets);
                 
            % 创建目标稀疏矩阵
            t_mat = sparse(double(target_rows), double(time_valid), double(tfr_vals), fLen, tLen);
            
            % 将稀疏矩阵转换为全矩阵并赋值给Ts
            Ts = gpuArray(full(t_mat));
        end
    end

    %% 归一化
    h0 = 1;
    Ts = Ts / N / h0;

    %% 返回结果
    Ts = gather(Ts);
    tfr = gather(tfr);
    omega = gather(omega);
end