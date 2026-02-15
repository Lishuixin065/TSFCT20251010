function [Tx, tfc0, t, f, T] = TET_2nd_gpu(x, fs, s)
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
    T = zeros(tLen, fLen, 'gpuArray');
    
    L0 = 1;
    gs = s;
    
    %% 向量化窗口计算和STFT
    % 预计算时间差矩阵 (N × N)
    t_diff_matrix = t' - tao;
    
    % 计算所有窗口函数的公共指数项
    exp_term = exp(-2*pi^2*gs^2*t_diff_matrix.^2./L0);
    
    % 预计算所有时间点的窗口函数
    gh0_all = 1./sqrt(L0) .* exp_term;
    gh1_all = -1i*2*pi*gs^2*t_diff_matrix.*(L0).^(-3/2).*exp_term;
    gh2_all = gs^2*((L0).^(-3/2)-((2*pi*gs*t_diff_matrix).^2).*(L0).^(-5/2)).*exp_term;
    
    %% 并行计算所有时间点的STFT
    for tidx = 1:N
        % 获取当前时间点的窗口
        gh0 = conj(gh0_all(tidx, :));
        gh1 = conj(gh1_all(tidx, :));
        gh2 = conj(gh2_all(tidx, :));
        
        % 计算FFT
        tf0 = fft(gh0 .* x);
        tf1 = fft(gh1 .* x);
        tf2 = fft(gh2 .* x);
        
        % 存储结果
        tfc0(tidx, :) = tf0(1:L);
        tfc1(tidx, :) = tf1(1:L);
        tfc2(tidx, :) = tf2(1:L);
    end
    
    %% 向量化计算群延迟
    M0 = tfc2 .* tfc0 - tfc1 .* tfc1;
    M1 = -tfc0 .* tfc1;
    
    % 创建时间索引矩阵 - 必须在GPU上
    time_indices = gpuArray(repmat((0:N-1)', 1, L));  % 已经是GPU
    Groupdelay = time_indices * del_t + M1 ./ M0 ./ (2.0*pi*1i);
    
    % 处理小值
    gamma = eps;
    small_mask = abs(tfc0) < gamma & abs(M0) < gamma;
    tfc0(small_mask) = 0;
    Groupdelay(small_mask) = 0;
    
    % 计算目标时间索引
    T = round(real(Groupdelay) ./ del_t);
    
    %% 优化的重分配 - 只保留m==b的情况（完全GPU实现）
    Tx = zeros(tLen, fLen, 'gpuArray');
    
    % 创建m==b的掩码（必须在GPU上创建索引矩阵）
    source_indices = repmat(gpuArray((1:tLen)'), 1, fLen);  % 关键修复！
    mask = (T == source_indices);
    
    % 直接赋值（由于只保留m==b，不需要累加）
    Tx(mask) = tfc0(mask);
    
    %% 将结果传回CPU
    Tx = gather(Tx);
    tfc0 = gather(tfc0);
    t = gather(t);
    f = gather(f);
    T = gather(T);
end