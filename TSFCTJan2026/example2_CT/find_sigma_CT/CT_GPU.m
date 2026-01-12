
function [ct,  f, c] = CT_gpu(x, fs, sigma, crange)
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

    % Pre-compute frequency base matrix
   
    % Pre-allocate outputs
    ct = gpuArray.zeros(nchirp, nfreq, N, 'single');
  
    % Pre-allocate workspace
    tf0_full = complex(gpuArray.zeros(nchirp, N, 'single'));
  
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
      
            time_sq_vals = time_sq(tau_idx);
            gauss_vals = gaussian_part(tau_idx);
            chirp_phase = exp(-1i * single(pi) * (c_col * time_sq_vals));
            g_window = gauss_vals .* chirp_phase;

            % Compute TF0, TF1, TF2
            x_segment_expanded = repmat(single(x_segment), nchirp, 1);
            tf0_seg = x_segment_expanded .* g_window;
       

            % Reset and fill
            tf0_full(:, circ_idx) = 0;
      

            tf0_full(:, circ_idx) = tf0_seg;
         

            % FFT
            tf0_fft = fft(tf0_full, [], 2) / single(fs);
       

            tf0 = tf0_fft(:, 1:nfreq);
     

            % Compute moments
           
            non_zero = abs(tf0) >= single(1e-10);

            ct_slice = gpuArray.zeros(nchirp, nfreq, 'single');
          
            ct_slice(non_zero) = tf0(non_zero);


            % Store results
            ct(:, :, tidx) = ct_slice;    
        end
    end

    % Final update
    fprintf('Processing iteration: %d/%d (100.0%%)\n', N, N);

    % Return results to CPU
    ct = gather(ct);
    f = gather(f);
    c = gather(c);
    fprintf('Extreme performance GPU computation completed!\n');
end



