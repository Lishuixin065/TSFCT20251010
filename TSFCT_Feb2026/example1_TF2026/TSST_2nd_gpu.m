function [Tx, tfc0, t, f, T] = TSST_2nd_gpu(x, fs, s)
    % 将输入数据转移到GPU
    x = gpuArray(x);
    
    [xrow, ~] = size(x);
    if (xrow ~= 1)
        x = x';
    end
    
    %% 预处理
    N = length(x);
    tao = gpuArray((0:N-1)/fs);
    t = gpuArray((0:N-1)/fs);
    
    if mod(N, 2) == 0
        f = gpuArray((0:N/2)*fs/N);
    else
        f = gpuArray((0:(N-1)/2)*fs/N);
    end
    
    L = length(f);
    del_t = 1/fs;
    tLen = N;
    fLen = L;
    
    %% 在GPU上预分配数组
    tfc0 = zeros(tLen, fLen, 'gpuArray');
    tfc1 = zeros(tLen, fLen, 'gpuArray');
    tfc2 = zeros(tLen, fLen, 'gpuArray');
 
    
    L0 = 1;
    gs = s;
    
    %% 直接在GPU上计算所有时间点的窗口函数
    % 创建时间差矩阵 (N × N)
    t_matrix = bsxfun(@minus, t(:), tao(:)');
    
    % 在GPU上计算所有窗口函数
    exp_term = exp(-2*pi^2*gs^2*t_matrix.^2./L0);
    
    gh0_all = 1./sqrt(L0) .* exp_term;
    gh1_all = -1i*2*pi*gs^2*t_matrix.*(L0).^(-3/2).*exp_term;
    gh2_all = gs^2*((L0).^(-3/2)-((2*pi*gs*t_matrix).^2).*(L0).^(-5/2)).*exp_term;
    
    %% 对每个时间点进行FFT计算（可以进一步向量化）
    for tidx = 1:N
        gh0 = conj(gh0_all(tidx, :));
        gh1 = conj(gh1_all(tidx, :));
        gh2 = conj(gh2_all(tidx, :));
        
        tf0 = fft(gh0 .* x);
        tf1 = fft(gh1 .* x);
        tf2 = fft(gh2 .* x);
        
        tfc0(tidx, :) = tf0(1:L);
        tfc1(tidx, :) = tf1(1:L);
        tfc2(tidx, :) = tf2(1:L);
    end
    
    %% 向量化计算群延迟
    M0 = tfc2 .* tfc0 - tfc1 .* tfc1;
    M1 = -tfc0 .* tfc1;
    
    time_indices = gpuArray(repmat((0:N-1)', 1, L));
    Groupdelay = time_indices * del_t + M1 ./ M0 ./ (2.0*pi*1i);
    
    % 处理小值
    gamma = eps;
    small_mask = abs(tfc0) < gamma & abs(M0) < gamma;
    tfc0(small_mask) = 0;
    Groupdelay(small_mask) = 0;
    
    T = round(real(Groupdelay) ./ del_t);
    
    %% ===== 优化的重分配 - 纯 GPU 向量化实现（保持原变量名）=====
    Tx = zeros(tLen, fLen, 'gpuArray');
    
    % 创建有效掩码
    valid_mask = T >= 1 & T <= tLen;
    
    if any(valid_mask(:))
        % 获取所有有效位置的索引
        [rows, cols] = find(valid_mask);  % rows = 目标时间索引，cols = 频率索引
        
        % 获取对应的值
        lin_idx = rows + (cols-1) * tLen;  % 原始线性索引
        tfc0_vals = tfc0(lin_idx);
        GroupDelay_vals = Groupdelay(lin_idx);  % 保持原变量名
        
        % 重新计算目标时间位置（从 GroupDelay 推导）
        target_rows = round(real(GroupDelay_vals) ./ del_t);  % 保持原逻辑
        
        % 确保目标时间位置在有效范围内
        valid_rows = target_rows >= 1 & target_rows <= tLen;  % 保持原变量名
        
        % 只保留有效的三元组
        cols = cols(valid_rows);
        tfc0_vals = tfc0_vals(valid_rows);
        target_rows = target_rows(valid_rows);
        
        % 构造新的线性索引（用于在 Tx 中的目标位置）
        target_lin_idx = target_rows + (cols - 1) * tLen;
        
        % 使用 groupsummary 在 GPU 上分组累加
        [unique_target_idx, ~, group_id] = unique(target_lin_idx);
        accum_vals = groupsummary(tfc0_vals, group_id, 'sum');  % 保持原逻辑名
        
        % 将累加结果填入 Tx
        Tx_flat = zeros(tLen * fLen, 1, 'gpuArray');
        Tx_flat(unique_target_idx) = accum_vals;
        Tx = reshape(Tx_flat, [tLen, fLen]);
    end
    
    %% 返回结果（转换到CPU）
    Tx = gather(Tx);
    tfc0 = gather(tfc0);
    t = gather(t);
    f = gather(f);
    T = gather(T);
end