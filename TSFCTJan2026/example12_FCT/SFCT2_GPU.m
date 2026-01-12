function [tfc0,  tfrtic, tcrtic, T, R] = SFCT2_GPU(x, fs, s, chrrange)
    %% GPU加速SFCT（无压缩版本）
    
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    
    N = length(x);
    tao = (0:N-1)/fs;   
    t = (0:N-1)/fs;
    
    %% === GPU数据传输 ===
    if gpuDeviceCount > 0
        x_gpu = gpuArray(single(x));
        t_gpu = gpuArray(single(t));
        tao_gpu = gpuArray(single(tao));
    else
        error('No GPU available!');
    end
    
    if mod(N,2)==0
        f = (0:N/2)*fs/N;
    else
        f = (0:(N-1)/2)*fs/N;
    end
    
    L = length(f);
    crate = linspace(-chrrange, chrrange, 2*round(N/2)+1);
    tfrtic = f;
    tcrtic = crate;
    cLen = length(crate);
    
    %% 预分配结果（留在GPU上）
    tfc0_gpu = gpuArray.zeros(cLen, N, L, 'single');
    tfc1_gpu = gpuArray.zeros(cLen, N, L, 'single');
    tfc2_gpu = gpuArray.zeros(cLen, N, L, 'single');
    T_gpu = gpuArray.zeros(cLen, N, L, 'single');
    R_gpu = gpuArray.zeros(cLen, N, L, 'single');
    
    del_t = single(1/fs);
    
    fprintf('GPU加速SFCT: N=%d, chirp rates=%d\n', N, cLen);
    tic_total = tic;
    
    %% 预计算tau矩阵
    tau_matrix = t_gpu' - tao_gpu;  % N×N
    
    %% === 主计算循环 ===
    for cidx = 1:cLen
        ahirp = single(crate(cidx));
        
        %% === 窗函数计算 ===
        L0 = 1 + 1i*2*pi*ahirp*s^2;
        t_sq = tau_matrix.^2;
        exp_arg = -2*pi^2*s^2*t_sq./L0;
        exp_part = exp(exp_arg);
        
        % 三种窗函数
        gh0_all = (1./sqrt(L0)) .* exp_part;
        gh1_all = -1i*2*pi*s^2*tau_matrix .* (L0).^(-3/2) .* exp_part;
        gh2_all = s^2 * ((L0).^(-3/2) - ((2*pi*s*tau_matrix).^2) .* (L0).^(-5/2)) .* exp_part;
        
        gh0_all = conj(gh0_all);
        gh1_all = conj(gh1_all);
        gh2_all = conj(gh2_all);
        
        %% === FFT计算 ===
        tf0_all = fft(bsxfun(@times, gh0_all, x_gpu), [], 2);
        tf1_all = fft(bsxfun(@times, gh1_all, x_gpu), [], 2);
        tf2_all = fft(bsxfun(@times, gh2_all, x_gpu), [], 2);
        
        tf0_all = tf0_all(:, 1:L);
        tf1_all = tf1_all(:, 1:L);
        tf2_all = tf2_all(:, 1:L);
        
        %% === 参数估计 ===
        M0 = tf2_all .* tf0_all - tf1_all .* tf1_all;
        M1 = -tf0_all .* tf1_all;
        M2 = tf0_all .* tf0_all;
        
        epsilon = single(1e-10);
        M0_abs = abs(M0);
        M0(M0_abs < epsilon) = epsilon;
        
        tidx_matrix = repmat((0:N-1)', 1, L);
        omega1 = tidx_matrix*del_t + real(M1 ./ M0 ./ (2.0*pi*1i));
        lambda1 = ahirp + real(M2 ./ M0 ./ (2.0*pi*1i));
        
        %% 存储结果（保持在GPU上）
        tfc0_gpu(cidx, :, :) = tf0_all;
        tfc1_gpu(cidx, :, :) = tf1_all;
        tfc2_gpu(cidx, :, :) = tf2_all;
        T_gpu(cidx, :, :) = omega1;
        R_gpu(cidx, :, :) = lambda1;
        
        % 显示进度
        if mod(cidx, max(1, floor(cLen/10))) == 0
            fprintf('进度: %d%%\n', round(cidx/cLen*100));
        end
    end
    
    %% === 传回CPU ===
    fprintf('传输结果回CPU...\n');
    tic_transfer = tic;
    
    tfc0 = gather(tfc0_gpu);
  
    T = gather(T_gpu);
    R = gather(R_gpu);
    
    fprintf('传输完成，耗时%.2f秒\n', toc(tic_transfer));
    
    total_time = toc(tic_total);
    fprintf('GPU加速SFCT完成！总时间: %.2f秒 (平均%.1f ms/chirp)\n', ...
            total_time, total_time/cLen*1000);
end