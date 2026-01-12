function [tfrsq1, ct, lam, omega, f, c] = CT_2nd_squeeze(x, fs, sigma, crange, num_squeeze)

    % 主函数：二阶时频变换GPU加速版（完整GPU版）
    % 输入：
    %   x - 输入信号
    %   fs - 采样频率
    %   sigma - 高斯窗标准差
    %   crange - 调频率范围
    %   num_squeeze - 压缩迭代次数（可选，默认不压缩）
    % 输出：
    %   ct, lam, omega, f, c - 原始CT变换结果
    %   tfrsq - 压缩后的时频表示（如果num_squeeze>0）
    
    % ========== 第一阶段：CT变换计算 ==========
    fprintf('========== CT二阶时频变换计算开始 ==========\n');
    
    % 设置默认压缩次数
    if nargin < 5
        num_squeeze = 0;
        do_squeeze = false;
    else
        do_squeeze = (num_squeeze > 0);
    end
    
    % GPU数据处理
    if ~isa(x, 'gpuArray')
        x = gpuArray(x);
    end
    if iscolumn(x), x = x'; end

    N = length(x);
    half_N = round(N/2);
    win_center = half_N + 1;

    % 预分配常量（单精度）
    c = gpuArray(single(linspace(-crange, crange, 2*half_N+1)));
    nchirp = length(c);
    c_col = c(:);

    % 窗口时间向量
    window_indices = -half_N:half_N;
    window_time = gpuArray(single(window_indices / fs));
    time_sq = single(window_time.^2);

    L = single(1/sigma/sqrt(2*pi));
    gaussian_part = single(L .* exp(-0.5*single(window_time.^2)/single(sigma^2)));

    % 频率向量
    if mod(N, 2) == 0
        f = single((0:N/2) * fs/N);
    else
        f = single((0:(N-1)/2) * fs/N);
    end
    f = gpuArray(f);
    nfreq = length(f);
    df = single(f(2) - f(1));

    % 预计算频率基矩阵
    freq_base_mat = gpuArray(single((0:(nfreq-1)) * df));
    freq_base_mat = repmat(freq_base_mat, [nchirp, 1]);

    % 预分配输出
    ct = gpuArray.zeros(nchirp, nfreq, N, 'single');
    lam = gpuArray.zeros(nchirp, nfreq, N, 'single');
    omega = gpuArray.zeros(nchirp, nfreq, N, 'single');

    % 预分配工作空间
    tf0_full = complex(gpuArray.zeros(nchirp, N, 'single'));
    tf1_full = complex(gpuArray.zeros(nchirp, N, 'single'));
    tf2_full = complex(gpuArray.zeros(nchirp, N, 'single'));

    fprintf('CT计算开始...\n');
    fprintf('信号长度: %d, 调频率数: %d, 频率数: %d\n', N, nchirp, nfreq);
    
    % 进度跟踪
    total_iterations = 0;
    display_interval = max(100, round(N/20));
    fprintf('进度: ');
    
    % 主计算循环
    for tidx = 1:N
        total_iterations = total_iterations + 1;
        
        % 进度显示
        if mod(total_iterations, display_interval) == 0 || tidx == N
            progress_percent = 100 * total_iterations / N;
            fprintf('%.1f%% ', progress_percent);
        end
        
        % 计算当前时间点的tau范围
        left_bound = -min([half_N-1, N, tidx-1]);
        right_bound = min([half_N-1, N, N-tidx]);
        tau = left_bound:right_bound;
        tau_len = length(tau);
        
        if tau_len == 0
            continue;
        end
        
        % 索引计算
        circ_idx = mod(N + tau, N) + 1;
        x_segment = x(tidx + tau);
        tau_idx = win_center + tau;
        
        % 计算窗口函数
        time_vals = window_time(tau_idx);
        time_sq_vals = time_sq(tau_idx);
        gauss_vals = gaussian_part(tau_idx);
        chirp_phase = exp(-1i * single(pi) * (c_col * time_sq_vals));
        g_window = gauss_vals .* chirp_phase;
        
        % 计算TF0, TF1, TF2
        x_segment_expanded = repmat(single(x_segment), nchirp, 1);
        tf0_seg = x_segment_expanded .* g_window;
        tf1_seg = tf0_seg .* time_vals;
        tf2_seg = tf0_seg .* time_sq_vals;
        
        % 重置并填充
        tf0_full(:, circ_idx) = 0;
        tf1_full(:, circ_idx) = 0;
        tf2_full(:, circ_idx) = 0;
        
        tf0_full(:, circ_idx) = tf0_seg;
        tf1_full(:, circ_idx) = tf1_seg;
        tf2_full(:, circ_idx) = tf2_seg;
        
        % FFT
        tf0_fft = fft(tf0_full, [], 2) / single(fs);
        tf1_fft = fft(tf1_full, [], 2) / single(fs);
        tf2_fft = fft(tf2_full, [], 2) / single(fs);
        
        tf0 = tf0_fft(:, 1:nfreq);
        tf1 = tf1_fft(:, 1:nfreq);
        tf2 = tf2_fft(:, 1:nfreq);
        
        % 计算矩
        M0 = tf0 .* tf2 - tf1 .* tf1;
        non_zero = abs(M0) >= single(1e-10);
        
        ct_slice = gpuArray.zeros(nchirp, nfreq, 'single');
        omega1 = gpuArray.zeros(nchirp, nfreq, 'single');
        lam1 = gpuArray.zeros(nchirp, nfreq, 'single');
        
        ct_slice(non_zero) = tf0(non_zero);
        
        if any(non_zero(:))
            M0_nz = M0(non_zero);
            M1_nz = tf0(non_zero) .* tf1(non_zero);
            M2_nz = -(tf0(non_zero) .* tf0(non_zero));
            
            inv_M0_nz = single(1) ./ M0_nz / single(2*pi);
            omega1(non_zero) = imag(M1_nz .* inv_M0_nz);
            lam1(non_zero) = imag(M2_nz .* inv_M0_nz);
        end
        
        % 存储结果
        ct(:, :, tidx) = ct_slice;
        lam(:, :, tidx) = c_col + lam1;
        omega(:, :, tidx) = omega1 + freq_base_mat;
    end
    
    fprintf('\nCT计算完成！\n');
    
    % ========== 第二阶段：可选的GPU压缩变换 ==========
    tfrsq1 = [];
    
    if do_squeeze
        fprintf('\n========== 开始GPU压缩变换 ==========\n');
        fprintf('压缩迭代次数: %d\n', num_squeeze);
        
        % 调用GPU压缩函数
        tfrsq1 = Msqueeze_CT_gpu_complete(x, ct, lam, omega, f, c, num_squeeze);
        
        fprintf('GPU压缩变换完成！\n');
    end
    
    % 传回CPU - 只在最后一次性传输所有数据
    ct = gather(ct);
    lam = gather(lam);
    omega = gather(omega);
    f = gather(f);
    c = gather(c);
    if do_squeeze
        tfrsq1 = gather(tfrsq1);
    end
    
    fprintf('========== 全部处理完成 ==========\n');
    
    % 清理GPU内存
    if gpuDeviceCount > 0
        reset(gpuDevice);
    end
end

% =========================================================================
% 附属函数：GPU压缩变换
% =========================================================================
function tfrsq1 =Msqueeze_CT_gpu_complete(x, ct, lam, omega, f, c, num_squeeze)
   
 
    % 确保数据在GPU上
    if ~isa(ct, 'gpuArray')
        ct = gpuArray(single(ct));
        lam = gpuArray(single(lam));
        omega = gpuArray(single(omega));
        f = gpuArray(single(f));
        c = gpuArray(single(c));
    end
    
    [Nc, Nf, Nt] = size(ct);
    
    % 计算分辨率
    df = f(2) - f(1);
    dc = c(2) - c(1);
    
    % 阈值
    if isa(x, 'gpuArray')
        x_cpu = gather(x);
    else
        x_cpu = x;
    end
    Threshold = single(0.0001 * mean(abs(x_cpu).^2));
    
    tfrsq_temp = ct;
    
    for iter = 1:num_squeeze
        fprintf('迭代 %d/%d: ', iter, num_squeeze);
        
        % 1. 一次性计算所有索引
        omega_idx = round(omega / df) + 1;
        lam_idx = round((lam + c(end)) / dc) + 1;
        
        % 2. 预分配输出
        tfrsq_out = gpuArray.zeros(Nc, Nf, Nt, 'single');
        
        % 3. 批量处理所有时间点
        progress_interval = max(1, round(Nt/20));
        
        for tt = 1:Nt
            if mod(tt, progress_interval) == 0
                fprintf('.');
            end
            
            % 提取当前时间片并展平
            ct_flat = reshape(tfrsq_temp(:, :, tt), [], 1);
            omega_flat = reshape(omega_idx(:, :, tt), [], 1);
            lam_flat = reshape(lam_idx(:, :, tt), [], 1);
            
            % 应用阈值
            valid_mask = abs(ct_flat) > Threshold;
            ct_valid = ct_flat(valid_mask);
            omega_valid = omega_flat(valid_mask);
            lam_valid = lam_flat(valid_mask);
            
            % 检查边界
            valid_idx = (omega_valid >= 1) & (omega_valid <= Nf) & ...
                        (lam_valid >= 1) & (lam_valid <= Nc);
            
            ct_valid = ct_valid(valid_idx);
            omega_valid = omega_valid(valid_idx);
            lam_valid = lam_valid(valid_idx);
            
            if isempty(ct_valid)
                continue;
            end
            
            % 方法1：使用histcounts2进行分箱（GPU版本）
            % 注意：histcounts2在GPU上可能不支持，使用替代方法
            
            % 创建线性索引
            linear_idx = (lam_valid - 1) * Nf + omega_valid;
            
            % 使用unique和手动累加（更高效的GPU实现）
            [unique_linear, ~, ic] = unique(linear_idx);
            
            % 预分配并累加
            accum = zeros(length(unique_linear), 1, 'single', 'gpuArray');
            
            % 批量累加（避免循环）
            for i = 1:length(unique_linear)
                mask = (ic == i);
                accum(i) = sum(ct_valid(mask));
            end
            
            % 重建矩阵
            tfrsq_slice = zeros(Nc, Nf, 'single', 'gpuArray');
            tfrsq_slice(unique_linear) = accum;
            tfrsq_out(:, :, tt) = tfrsq_slice;
        end
        
        fprintf('\n');
        tfrsq_temp = tfrsq_out;
    end
    
    tfrsq1 = tfrsq_temp;
end