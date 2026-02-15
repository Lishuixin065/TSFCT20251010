function [Ts, tfr0, tfrtic] = SST_2nd_hop_gpu(x, sigma, hop, Hz)
    % GPU加速的带hop的SST第二版
    % 注意参数顺序：x, sigma, hop, Hz

    [xrow, xcol] = size(x);
    if (xcol ~= 1)
        error('x must have only one column');
    end
    
    if gpuDeviceCount == 0
        warning('GPU not available, running on CPU');
        [Ts, tfr0, tfrtic] = SST_2nd_hop(x, sigma, hop, Hz);
        return;
    end
    
    % 重置GPU
    gpuDevice(1);
    
    N = xrow;
    t = 1:hop:N;
    L = N / Hz;  % 信号总时长
    
    % 频率向量
    if mod(N, 2) == 0
        f = (0:hop:round(N/2)) * Hz / N;
    else
        f = (0:hop:round((N-1)/2)) * Hz / N;
    end
    tfrtic = f;
    
    tLen = length(t);
    fLen = length(tfrtic);
    
    % 关键：time的计算（与原始代码一致）
    ht = -tLen:tLen;
    time = L * ht / tLen / hop;  % 注意这里的计算
    
    %% 数据转移到GPU
    x_gpu = gpuArray(single(x));
    time_gpu = gpuArray(single(time));
    
    %% 预计算窗函数（与原始代码一致）
    sigma_sq = single(sigma^2);
    factor = single(-2 * pi^2 * sigma_sq);
    
    time_sq = time_gpu.^2;
    exp_term = exp(factor * time_sq);
    
    % 注意：原始代码使用列向量，所以使用transpose
    h = exp_term';        % 高斯窗
    th = (time_gpu .* exp_term)';    % t * 高斯窗
    tth = (time_sq .* exp_term)';    % t^2 * 高斯窗
    
    %% STFT计算（GPU优化）
    tfr0_gpu = complex(gpuArray.zeros(fLen, tLen, 'single'));
    tfr1_gpu = tfr0_gpu;
    tfr2_gpu = tfr0_gpu;
    
    % 预分配工作数组
    tf0_work = complex(gpuArray.zeros(tLen, 1, 'single'));
    tf1_work = tf0_work;
    tf2_work = tf0_work;
    
    % 主循环
    for tidx = 1:tLen
        ti = t(tidx);
        
        % 计算tau范围（与原始代码一致）
        Lh = tLen;  % 注意：Lh = tLen
        tau = -min([round(tLen/2)-1, Lh, ti-1]):min([round(tLen/2)-1, Lh, N-ti]);
        
        % 转换为列向量
        tau = tau(:);
        
        % 计算indices
        indices = rem(tLen + tau, tLen) + 1;
        indices = indices(:);
        
        % 重置工作数组
        tf0_work(:) = 0;
        tf1_work(:) = 0;
        tf2_work(:) = 0;
        
        % 获取窗函数值（注意：Lh+1+tau）
        tau_idx = Lh + 1 + tau;
        h_vals = h(tau_idx);
        th_vals = th(tau_idx);
        tth_vals = tth(tau_idx);
        
        % 获取信号片段
        x_segment = x_gpu(ti + tau);
        
        % 确保所有向量都是列向量
        h_vals = h_vals(:);
        th_vals = th_vals(:);
        tth_vals = tth_vals(:);
        x_segment = x_segment(:);
        
        % 赋值到工作数组
        tf0_work(indices) = x_segment .* h_vals;
        tf1_work(indices) = x_segment .* th_vals;
        tf2_work(indices) = x_segment .* tth_vals;
        
        % FFT计算
        tf0_fft = fft(tf0_work);
        tf1_fft = fft(tf1_work);
        tf2_fft = fft(tf2_work);
        
        % 存储结果（只取前fLen个点）
        tfr0_gpu(:, tidx) = tf0_fft(1:fLen);
        tfr1_gpu(:, tidx) = tf1_fft(1:fLen);
        tfr2_gpu(:, tidx) = tf2_fft(1:fLen);
    end
    
    %% 频率重分配计算
    del_f = tfrtic(2) - tfrtic(1);
    del_f_gpu = gpuArray(single(del_f));
    
    % 避免除以零
    denominator = tfr1_gpu .* tfr1_gpu - tfr0_gpu .* tfr2_gpu;
    denominator(denominator == 0) = eps('single');
    
    % 计算瞬时频率
    omega1 = real((tfr0_gpu .* tfr1_gpu) ./ denominator ./ (2.0 * pi * 1i));
    omega1 = round(omega1 ./ del_f_gpu);
    
    % 计算omega（注意：+1补偿，与原代码一致）
    eta_vec = gpuArray(single(1:fLen)');  % 在GPU上创建
    omega_gpu = eta_vec - omega1 + 1;
    
    %% ===== 同步赋值（完全GPU优化）=====
    Ts_gpu = gpuArray.zeros(fLen, tLen, 'single');
    eps_val = single(0.0001);  % 与原代码阈值一致
    
    % 使用高效的向量化方法
    valid_mask = abs(tfr0_gpu) > eps_val;
    
    if any(valid_mask(:))
        [rows, cols] = find(valid_mask);
        
        % 使用线性索引
        lin_idx = rows + (cols-1) * fLen;
        
        tfr_vals = tfr0_gpu(lin_idx);
        omega_vals = omega_gpu(lin_idx);
        
        % 计算目标位置（取整）
        target_rows = round(omega_vals);
        
        % 筛选有效位置
        valid_target = (target_rows >= 1) & (target_rows <= fLen);
        
        if any(valid_target)
            valid_rows = target_rows(valid_target);
            valid_cols = cols(valid_target);
            valid_tfr_vals = tfr_vals(valid_target);
            
            % 创建目标线性索引
            target_lin_idx = valid_rows + (valid_cols-1) * fLen;
            
            % === 修正：使用 groupsummary 在 GPU 上累加 ===
            [unique_target_idx, ~, group_id] = unique(target_lin_idx);
            sum_vals = groupsummary(valid_tfr_vals, group_id, 'sum');
            
            % 将累加结果填入 Ts_gpu
            Ts_flat = gpuArray.zeros(fLen * tLen, 1, 'single');
            Ts_flat(unique_target_idx) = sum_vals;
            Ts_gpu = reshape(Ts_flat, [fLen, tLen]);
        end
    end
    
    %% 最终归一化（注意：除以N，不是tLen）
    Ts_gpu = Ts_gpu / N;
    
    %% 传回结果到CPU
    Ts = gather(Ts_gpu);
    tfr0 = gather(tfr0_gpu);
end