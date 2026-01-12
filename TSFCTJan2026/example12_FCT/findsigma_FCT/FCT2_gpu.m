function [tfc0] = FCT2_gpu(x,fs,s,gddRange)
% % Inputs:
%   x         : Input signal (vector)
%   fs        : Sampling frequency (Hz)
%   s         : Standard deviation of Gaussian window (in samples)
%   chrrange  : Maximum chirp rate range (determines GDD axis resolution)
%
% Outputs:
%   tfc0      : FCT coefficients (using base Gaussian window g)

    
    %% === 1. Data Preparation ===
    [xrow, ~] = size(x);
    if xrow ~= 1
        x = x';
    end
    
    N = length(x);
    tau = (0:N-1)/fs;   
    t = (0:N-1)/fs;
    
    if gpuDeviceCount == 0
        error('No GPU available!');
    end
    
    % Frequency axis
    if mod(N, 2) == 0
        f = (0:N/2) * fs / N;
    else
        f = (0:(N-1)/2) * fs / N;
    end
    
    Lf = length(f);
    gddValues = linspace(-gddRange, gddRange, 2*round(N/2)+1);

    Ngdd = length(gddValues);
    


    
    %% === 2. Single GPU Transfer ===
    fprintf('Transferring all data to GPU...\n');
    x_gpu = gpuArray(single(x));
    t_gpu = gpuArray(single(t));
    tau_gpu = gpuArray(single(tau));
    
    % Precompute time difference matrix
    dtMat = t_gpu' - tau_gpu;
    
    %% === 3. SFFT Computation (all data remains on GPU) ===
    fprintf('Starting SFFT computation...\n');
    ticSFFT = tic;
    
    % Preallocate GPU arrays
    tf0_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
   
    

    
    for gddIdx = 1:Ngdd
        gddVal = single(gddValues(gddIdx));
        
        % Window function computation with GDD parameter
        L0 = 1 + 1i*2*pi*gddVal*s^2;
        t2 = dtMat.^2;
        expArg = -2*pi^2*s^2*t2./L0;
        expTerm = exp(expArg);
        
        % Window and its derivatives for GDD estimation
        g0 = (1./sqrt(L0)) .* expTerm;
     
        % Complex conjugate
        g0 = conj(g0);
      
        % Fourier transform
        tf0 = fft(bsxfun(@times, g0, x_gpu), [], 2);
     
        % Keep positive frequencies
        tf0 = tf0(:, 1:Lf);
      
        % Parameter estimation for time and GDD
      
     
       
        % Store results
        tf0_gpu(gddIdx, :, :) = tf0;
      
    end
     tfc0 = gather(tf0_gpu);
    fprintf('SFFT computation completed: %.2f seconds\n', toc(ticSFFT));
    
  
end