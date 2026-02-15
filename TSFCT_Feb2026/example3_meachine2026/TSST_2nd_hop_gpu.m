function [Ts, tfr0, t0, f, GroupDelay] = TSST_2nd_hop_gpu(x, Hz, sigma, hop)
    % GPU加速的TSST_2nd_hop版本 - 无进度显示
    % 包含相位补偿和群延迟计算

    [xrow, xcol] = size(x);
    if (xcol ~= 1)
        error('x must have only one column');
    end
    
    if gpuDeviceCount == 0
        [Ts, tfr0, tfr1, tfr2, t0, f, GroupDelay] = TSST_2nd_hop(x, Hz, sigma, hop);
        return;
    end
    
    % 重置GPU
    gpuDevice(1);
    
    N = xrow;
    t = 1:hop:N;
    L = N / Hz;
    
    % 频率向量
    if mod(N, 2) == 0
        f = (0:hop:round(N/2)) * Hz / N;
    else
        f = (0:hop:round((N-1)/2)) * Hz / N;
    end
    
    tLen = length(t);
    fLen = length(f);
    
    % 时间向量
    ht = single(-tLen:tLen);
    time = L * ht / tLen / hop;
    
    %% 数据转移到GPU
    x_gpu = gpuArray(single(x));
    time_gpu = gpuArray(single(time));
    f_gpu = gpuArray(single(f));
    
    %% 预计算窗函数
    sigma_sq = single(sigma^2);
    factor = single(-2 * pi^2 * sigma_sq);
    
    time_sq = time_gpu.^2;
    exp_term = exp(factor * time_sq);
    
    % 转置以匹配原始代码
    h = exp_term';
    th = (time_gpu .* exp_term)';
    tth = (time_sq .* exp_term)';
    
    %% STFT计算
    tfr0_gpu = complex(gpuArray.zeros(fLen, tLen, 'single'));
    tfr1_gpu = tfr0_gpu;
    tfr2_gpu = tfr0_gpu;
    
    % 工作数组
    tf0_work = complex(gpuArray.zeros(tLen, 1, 'single'));
    tf1_work = tf0_work;
    tf2_work = tf0_work;
    
    Lh = tLen;
    max_offset = round(tLen/2) - 1;
    
    for tidx = 1:tLen
        ti = t(tidx);
        
        % 计算tau范围
        left = min([max_offset, Lh, ti-1]);
        right = min([max_offset, Lh, N-ti]);
        tau = (-left:right)';
        
        % 计算索引
        indices = mod(tLen + tau, tLen) + 1;
        
        % 重置工作数组
        tf0_work(:) = 0;
        tf1_work(:) = 0;
        tf2_work(:) = 0;
        
        % 获取窗函数值
        tau_idx = Lh + 1 + tau;
        h_vals = h(tau_idx);
        th_vals = th(tau_idx);
        tth_vals = tth(tau_idx);
        x_vals = x_gpu(ti + tau);
        
        % 赋值
        tf0_work(indices) = x_vals .* h_vals;
        tf1_work(indices) = x_vals .* th_vals;
        tf2_work(indices) = x_vals .* tth_vals;
        
        % FFT
        tf0_fft = fft(tf0_work);
        tf1_fft = fft(tf1_work);
        tf2_fft = fft(tf2_work);
        
        % 存储
        tfr0_gpu(:, tidx) = tf0_fft(1:fLen);
        tfr1_gpu(:, tidx) = tf1_fft(1:fLen);
        tfr2_gpu(:, tidx) = tf2_fft(1:fLen);
    end
    
    %% 相位补偿
    % 时间向量t0
    t0_gpu = gpuArray(single(t / Hz));
    
    % 预计算相位因子
    phase_factor = exp(-2 * pi * 1i * f_gpu(:) * t0_gpu(:)');
    
    % 相位补偿（向量化）
    tfr0_gpu = tfr0_gpu .* phase_factor;
    tfr1_gpu = tfr1_gpu .* phase_factor;
    tfr2_gpu = tfr2_gpu .* phase_factor;
    
    %% 计算tfr1和tfr2的修正
    % 预计算常数
    sigma2_inv = 1 / sigma_sq;
    sigma4_inv = sigma2_inv^2;
    factor1 = -1 / (2 * pi * 1i);
    factor2 = factor1^2;
    
    % 计算修正（向量化）
    tfr1_gpu = -tfr1_gpu .* sigma2_inv .* factor1;
    tfr2_gpu = (-tfr0_gpu .* sigma2_inv + tfr2_gpu .* sigma4_inv) .* factor2;
    
    %% 计算群延迟（GroupDelay）
    % 计算分母和分子
    Denominator = tfr2_gpu .* tfr0_gpu - tfr1_gpu .* tfr1_gpu;
    Denominator(Denominator == 0) = eps('single');
    
    Numerator = tfr0_gpu .* tfr1_gpu;
    
    % 计算p
    p = Numerator ./ Denominator .* single(factor1);
    
    % 添加时间偏移（向量化）
    del_t = hop / Hz;
    tau_vec = gpuArray(single(1:tLen));
    p = bsxfun(@plus, p, tau_vec * del_t);
    
    % 计算群延迟并取整
    GroupDelay_gpu = round(real(p) / del_t);
    
    %% ===== 时间同步压缩（完全GPU实现）=====
    Ts_gpu = gpuArray.zeros(fLen, tLen, 'single');
    
    % 找到所有有效位置
    valid_mask = (GroupDelay_gpu >= 1) & (GroupDelay_gpu <= tLen);
    
    if any(valid_mask(:))
        % 获取所有有效位置的索引
        [rows, cols] = find(valid_mask);
        
        % 使用线性索引
        lin_idx = rows + (cols-1) * fLen;
        
        % 获取对应的值
        tfr0_vals = tfr0_gpu(lin_idx);
        GroupDelay_vals = GroupDelay_gpu(lin_idx);
        
        % 计算目标列位置（取整）
        target_cols = round(GroupDelay_vals);
        
        % 创建目标线性索引
        target_lin_idx = rows + (target_cols-1) * fLen;
        
        % === 修正：使用 groupsummary 在 GPU 上累加 ===
        [unique_target_idx, ~, group_id] = unique(target_lin_idx);
        sum_vals = groupsummary(tfr0_vals, group_id, 'sum');
        
        % 将累加结果填入 Ts_gpu
        Ts_flat = gpuArray.zeros(fLen * tLen, 1, 'single');
        Ts_flat(unique_target_idx) = sum_vals;
        Ts_gpu = reshape(Ts_flat, [fLen, tLen]);
    end
    
    %% 将结果传回CPU
    Ts = gather(Ts_gpu);
    tfr0 = gather(tfr0_gpu);
    t0 = gather(t0_gpu);
    f = gather(f_gpu);
    GroupDelay = gather(GroupDelay_gpu);
end