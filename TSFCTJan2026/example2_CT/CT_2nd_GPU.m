
function [ct, lambda, omega, f, c] = CT_2nd_GPU(x, fs, sigma, crange)
    % Final optimized version - GPU accelerated with single precision

    if ~isa(x, 'gpuArray')
        x = gpuArray(x);
    end
    if iscolumn(x), x = x'; end

    N = length(x);
    half_N = round(N/2);
    win_center = half_N + 1;

    % Pre-allocate constants with single precision
    c = gpuArray(single(linspace(-crange, crange, 2*half_N+1)));
    nchirp = length(c);
    c_col = c(:);

    % Window time vector
    window_indices = -half_N:half_N;
    window_time = gpuArray(single(window_indices / fs));
    time_sq = single(window_time.^2);

    L = single(1/sigma/sqrt(2*pi));
    gaussian_part = single(L .* exp(-0.5*single(window_time.^2)/single(sigma^2)));

    % Frequency vector
    if mod(N, 2) == 0
        f = single((0:N/2) * fs/N);
    else
        f = single((0:(N-1)/2) * fs/N);
    end
    f = gpuArray(f);
    nfreq = length(f);
    df = single(f(2) - f(1));

    % Pre-compute frequency base matrix
    freq_base_mat = gpuArray(single((0:(nfreq-1)) * df));
    freq_base_mat = repmat(freq_base_mat, [nchirp, 1]);

    % Pre-allocate outputs
    ct = gpuArray.zeros(nchirp, nfreq, N, 'single');
    lambda = gpuArray.zeros(nchirp, nfreq, N, 'single');
    omega = gpuArray.zeros(nchirp, nfreq, N, 'single');

    % Pre-allocate workspace
    tf0_full = complex(gpuArray.zeros(nchirp, N, 'single'));
    tf1_full = complex(gpuArray.zeros(nchirp, N, 'single'));
    tf2_full = complex(gpuArray.zeros(nchirp, N, 'single'));

    fprintf('Starting extreme performance GPU computation...\n');

    % Process in chunks for better memory management
    chunk_size = min(500, N);
    num_chunks = ceil(N / chunk_size);

    total_iterations = 0;  % 记录总处理次数
    display_interval = 100;  % 每100次显示一次

    for chunk = 1:num_chunks
        start_tidx = (chunk-1) * chunk_size + 1;
        end_tidx = min(chunk * chunk_size, N);

        for tidx = start_tidx:end_tidx
            total_iterations = total_iterations + 1;

            % Show progress every 100 iterations
            if mod(total_iterations, display_interval) == 1
                progress_percent = 100 * (total_iterations-1) / N;
                fprintf('Processing iteration: %d/%d (%.1f%%)\n', total_iterations-1, N, progress_percent);
            end

            % Calculate tau range for current time point
            left_bound = -min([half_N-1, N, tidx-1]);
            right_bound = min([half_N-1, N, N-tidx]);
            tau = left_bound:right_bound;
            tau_len = length(tau);

            if tau_len == 0
                continue;
            end

            circ_idx = mod(N + tau, N) + 1;
            x_segment = x(tidx + tau);
            tau_idx = win_center + tau;

            % Compute with single precision
            time_vals = window_time(tau_idx);
            time_sq_vals = time_sq(tau_idx);
            gauss_vals = gaussian_part(tau_idx);
            chirp_phase = exp(-1i * single(pi) * (c_col * time_sq_vals));
            g_window = gauss_vals .* chirp_phase;

            % Compute TF0, TF1, TF2
            x_segment_expanded = repmat(single(x_segment), nchirp, 1);
            tf0_seg = x_segment_expanded .* g_window;
            tf1_seg = tf0_seg .* time_vals;
            tf2_seg = tf0_seg .* time_sq_vals;

            % Reset and fill
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

            % Compute moments
            M0 = tf0 .* tf2 - tf1 .* tf1;
            non_zero = abs(M0) >= single(1e-10);

            ct_slice = gpuArray.zeros(nchirp, nfreq, 'single');
            omega1 = gpuArray.zeros(nchirp, nfreq, 'single');
            lambda1 = gpuArray.zeros(nchirp, nfreq, 'single');

            ct_slice(non_zero) = tf0(non_zero);

            if any(non_zero(:))
                M0_nz = M0(non_zero);
                M1_nz = tf0(non_zero) .* tf1(non_zero);
                M2_nz = -(tf0(non_zero) .* tf0(non_zero));

                inv_M0_nz = single(1) ./ M0_nz / single(2*pi);
                omega1(non_zero) = imag(M1_nz .* inv_M0_nz);
                lambda1(non_zero) = imag(M2_nz .* inv_M0_nz);
            end

            % Store results
            ct(:, :, tidx) = ct_slice;
            lambda(:, :, tidx) = c_col + lambda1;
            omega(:, :, tidx) = omega1 + freq_base_mat;
        end
    end

    % Final update
    fprintf('Processing iteration: %d/%d (100.0%%)\n', N, N);

    % Return results to CPU
    ct = gather(ct);
    lambda = gather(lambda);
    omega = gather(omega);
    f = gather(f);
    c = gather(c);
    fprintf('Extreme performance GPU computation completed!\n');
end



%% 
% function [ct, lambda, omega, f, c] = CT_2nd_GPU_optimized(x, fs, sigma, crange)
%     % 终极优化版本 - GPU加速 + 并行计算 + 内存优化
% 
%     % ============== 1. 初始化和数据准备 ==============
%     tic_init = tic;
% 
%     if ~isa(x, 'gpuArray')
%         x = gpuArray(single(x));
%     else
%         x = single(x);
%     end
% 
%     if iscolumn(x), x = x'; end
% 
%     N = length(x);
%     half_N = round(N/2);
%     win_center = half_N + 1;
% 
%     % 预计算常数
%     c = gpuArray(single(linspace(-crange, crange, 2*half_N+1)));
%     nchirp = length(c);
%     c_col = c(:);
% 
%     % 窗口时间向量
%     window_indices = -half_N:half_N;
%     window_time = gpuArray(single(window_indices / fs));
%     time_sq = window_time.^2;
% 
%     % 高斯窗口预计算
%     sigma_sq = sigma^2;
%     gaussian_part = single(1/(sigma*sqrt(2*pi))) .* exp(-0.5 * window_time.^2 / sigma_sq);
%     gaussian_part = gpuArray(gaussian_part);
% 
%     % 频率向量
%     if mod(N, 2) == 0
%         f = gpuArray(single((0:N/2) * fs/N));
%     else
%         f = gpuArray(single((0:(N-1)/2) * fs/N));
%     end
%     nfreq = length(f);
%     df = f(2) - f(1);
% 
%     % 预计算频率基础矩阵
%     freq_base = gpuArray(single(0:(nfreq-1))) * df;
%     freq_base_mat = repmat(freq_base, [nchirp, 1]);
% 
%     % 预分配输出
%     ct = gpuArray.zeros(nchirp, nfreq, N, 'single');
%     lambda = gpuArray.zeros(nchirp, nfreq, N, 'single');
%     omega = gpuArray.zeros(nchirp, nfreq, N, 'single');
% 
%     fprintf('初始化完成，耗时: %.2f秒\n', toc(tic_init));
% 
%     % ============== 2. 预计算循环不变量 ==============
%     tic_precomp = tic;
% 
%     % 预计算所有可能的tau索引和对应的circ_idx
%     tau_all = cell(N, 1);
%     circ_idx_all = cell(N, 1);
%     win_tau_idx_all = cell(N, 1);
% 
%     for tidx = 1:N
%         left_bound = -min([half_N-1, N, tidx-1]);
%         right_bound = min([half_N-1, N, N-tidx]);
%         tau = left_bound:right_bound;
%         tau_all{tidx} = tau;
%         circ_idx_all{tidx} = mod(N + tau, N) + 1;
%         win_tau_idx_all{tidx} = win_center + tau;
%     end
% 
%     % 预计算窗函数值矩阵（避免重复计算）
%     window_time_mat = repmat(window_time', 1, nchirp);
%     time_sq_mat = repmat(time_sq', 1, nchirp);
%     gaussian_mat = repmat(gaussian_part', 1, nchirp);
% 
%     fprintf('预计算完成，耗时: %.2f秒\n', toc(tic_precomp));
% 
%     % ============== 3. 主计算循环（优化版本） ==============
%     tic_main = tic;
% 
%     % 使用更大的chunk_size减少循环开销
%     chunk_size = min(1000, N);
%     num_chunks = ceil(N / chunk_size);
% 
%     % 预分配工作空间（在chunk级别）
%     tf0_full = complex(gpuArray.zeros(nchirp, N, 'single'));
%     tf1_full = complex(gpuArray.zeros(nchirp, N, 'single'));
%     tf2_full = complex(gpuArray.zeros(nchirp, N, 'single'));
% 
%     fprintf('开始主计算...\n');
% 
%     for chunk = 1:num_chunks
%         chunk_start = (chunk-1) * chunk_size + 1;
%         chunk_end = min(chunk * chunk_size, N);
%         chunk_len = chunk_end - chunk_start + 1;
% 
%         % 批量处理当前chunk的所有时间点
%         for t_local = 1:chunk_len
%             tidx = chunk_start + t_local - 1;
% 
%             % 进度显示
%             if mod(tidx, 500) == 0
%                 fprintf('处理进度: %d/%d (%.1f%%)\n', tidx, N, 100*tidx/N);
%             end
% 
%             % 获取预计算的索引
%             tau = tau_all{tidx};
%             if isempty(tau)
%                 continue;
%             end
% 
%             circ_idx = circ_idx_all{tidx};
%             win_tau_idx = win_tau_idx_all{tidx};
%             tau_len = length(tau);
% 
%             % 获取信号片段
%             x_segment = x(tidx + tau);
% 
%             % ============== 向量化计算核心部分 ==============
%             % 获取预计算的窗口值
%             time_vals = window_time(win_tau_idx);
%             time_sq_vals = time_sq(win_tau_idx);
%             gauss_vals = gaussian_part(win_tau_idx);
% 
%             % 计算chirp相位（向量化）
%             % 使用bsxfun进行矩阵乘法，避免显式循环
%             chirp_phase = exp(-1i * pi * (c_col * time_sq_vals));
% 
%             % 计算窗函数
%             g_window = gauss_vals .* chirp_phase;
% 
%             % 扩展信号片段
%             x_expanded = repmat(x_segment, nchirp, 1);
% 
%             % 计算TF0, TF1, TF2
%             tf0_seg = x_expanded .* g_window;
%             tf1_seg = tf0_seg .* time_vals;
%             tf2_seg = tf0_seg .* time_sq_vals;
% 
%             % 重置工作空间（只重置需要的部分）
%             tf0_full(:, circ_idx) = 0;
%             tf1_full(:, circ_idx) = 0;
%             tf2_full(:, circ_idx) = 0;
% 
%             tf0_full(:, circ_idx) = tf0_seg;
%             tf1_full(:, circ_idx) = tf1_seg;
%             tf2_full(:, circ_idx) = tf2_seg;
% 
%             % FFT计算（批量）
%             tf0_fft = fft(tf0_full, [], 2) / fs;
%             tf1_fft = fft(tf1_full, [], 2) / fs;
%             tf2_fft = fft(tf2_full, [], 2) / fs;
% 
%             tf0 = tf0_fft(:, 1:nfreq);
%             tf1 = tf1_fft(:, 1:nfreq);
%             tf2 = tf2_fft(:, 1:nfreq);
% 
%             % ============== 快速计算矩量 ==============
%             % 计算M0, M1, M2
%             M0 = tf0 .* tf2 - tf1 .* tf1;
%             M1 = tf0 .* tf1;
%             M2 = -tf0 .* tf0;
% 
%             % 使用逻辑索引避免条件判断
%             valid_mask = abs(M0) >= 1e-10;
% 
%             % 存储结果
%             ct_slice = zeros(nchirp, nfreq, 'single', 'gpuArray');
%             lambda_slice = zeros(nchirp, nfreq, 'single', 'gpuArray');
%             omega_slice = zeros(nchirp, nfreq, 'single', 'gpuArray');
% 
%             ct_slice(valid_mask) = tf0(valid_mask);
% 
%             if any(valid_mask(:))
%                 M0_valid = M0(valid_mask);
%                 M1_valid = M1(valid_mask);
%                 M2_valid = M2(valid_mask);
% 
%                 inv_M0 = 1 ./ (2 * pi * M0_valid);
% 
%                 omega_slice(valid_mask) = imag(M1_valid .* inv_M0);
%                 lambda_slice(valid_mask) = imag(M2_valid .* inv_M0);
%             end
% 
%             % 存储到最终数组
%             ct(:, :, tidx) = ct_slice;
%             lambda(:, :, tidx) = c_col + lambda_slice;
%             omega(:, :, tidx) = omega_slice + freq_base_mat;
%         end
%     end
% 
%     fprintf('主计算完成，耗时: %.2f秒\n', toc(tic_main));
% 
%     % ============== 4. 结果返回 ==============
%     tic_gather = tic;
% 
%     ct = gather(ct);
%     lambda = gather(lambda);
%     omega = gather(omega);
%     f = gather(f);
%     c = gather(c);
% 
%     fprintf('数据回传完成，耗时: %.2f秒\n', toc(tic_gather));
%     fprintf('总计算时间: %.2f秒\n', toc(tic_init));
% end