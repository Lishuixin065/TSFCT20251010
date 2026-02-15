function [tfComp, tfOrig, fAxis, gddAxis] = SFCT2_with_squeeze_gpu(x, fs, s, gddRange, sqNum, sqThr)
    % SFCT2_WITH_SQUEEZE_GPU - Sparse Fractional Fourier Transform with squeezing on GPU
    %
    % INPUTS:
    %   x       - Input signal (1D vector)
    %   fs      - Sampling frequency (Hz)
    %   s       - Window parameter (sigma)
    %   gddRange- Group delay dispersion range (maximum absolute GDD)
    %   sqNum   - Number of squeezing iterations
    %   sqThr   - Squeezing threshold factor
    %
    % OUTPUTS:
    %   tfComp  - Compressed time-frequency representation
    %   tfOrig  - Original time-frequency representation
    %   fAxis   - Frequency axis (Hz)
    %   gddAxis - Group delay dispersion axis (s^2)
    %
    % DESCRIPTION:
    %   This function computes the Sparse Fractional Fourier Transform (FCT)
    %   with iterative squeezing compression entirely on GPU. The transform
    %   considers group delay dispersion (GDD) as the second-order parameter.
    %
    % AUTHOR:
    %   [Your Name/Institution]
    %
    % VERSION:
    %   1.0 - [Date]
    %
    % REFERENCES:
    %   [1] Reference 1
    %   [2] Reference 2
    
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
    fAxis = f;
    gddAxis = gddValues;
    Ngdd = length(gddValues);
    
    fprintf('Integrated GPU pipeline: N=%d, GDD points=%d, squeeze iterations=%d\n', ...
            N, Ngdd, sqNum);
    ticTotal = tic;
    
    %% === 2. Single GPU Transfer ===
    fprintf('Transferring all data to GPU...\n');
    x_gpu = gpuArray(single(x));
    t_gpu = gpuArray(single(t));
    tau_gpu = gpuArray(single(tau));
    
    % Precompute time difference matrix
    dtMat = t_gpu' - tau_gpu;
    
    %% === 3. FCT Computation (all data remains on GPU) ===
    fprintf('Starting FCT computation...\n');
    ticFCT = tic;
    
    % Preallocate GPU arrays
    tf0_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
    tf1_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
    tf2_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
    timeEst_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
    gddEst_gpu = gpuArray.zeros(Ngdd, N, Lf, 'single');
    
    dt = single(1/fs);
    
    for gddIdx = 1:Ngdd
        gddVal = single(gddValues(gddIdx));
        
        % Window function computation with GDD parameter
        L0 = 1 + 1i*2*pi*gddVal*s^2;
        t2 = dtMat.^2;
        expArg = -2*pi^2*s^2*t2./L0;
        expTerm = exp(expArg);
        
        % Window and its derivatives for GDD estimation
        g0 = (1./sqrt(L0)) .* expTerm;
        g1 = -1i*2*pi*s^2*dtMat .* (L0).^(-3/2) .* expTerm;
        g2 = s^2 * ((L0).^(-3/2) - ((2*pi*s*dtMat).^2) .* (L0).^(-5/2)) .* expTerm;
        
        % Complex conjugate
        g0 = conj(g0);
        g1 = conj(g1);
        g2 = conj(g2);
        
        % Fourier transform
        tf0 = fft(bsxfun(@times, g0, x_gpu), [], 2);
        tf1 = fft(bsxfun(@times, g1, x_gpu), [], 2);
        tf2 = fft(bsxfun(@times, g2, x_gpu), [], 2);
        
        % Keep positive frequencies
        tf0 = tf0(:, 1:Lf);
        tf1 = tf1(:, 1:Lf);
        tf2 = tf2(:, 1:Lf);
        
        % Parameter estimation for time and GDD
        M0 = tf2 .* tf0 - tf1 .* tf1;
        M1 = -tf0 .* tf1;
        M2 = tf0 .* tf0;
        
        % Avoid division by zero
        epsVal = single(1e-10);
        M0Abs = abs(M0);
        M0(M0Abs < epsVal) = epsVal;
        
        % Time and GDD estimation
        tIdxMat = repmat((0:N-1)', 1, Lf);
        timeEst = tIdxMat*dt + real(M1 ./ M0 ./ (2.0*pi*1i));
        gddEst = gddVal + real(M2 ./ M0 ./ (2.0*pi*1i));
        
        % Store results
        tf0_gpu(gddIdx, :, :) = tf0;
        tf1_gpu(gddIdx, :, :) = tf1;
        tf2_gpu(gddIdx, :, :) = tf2;
        timeEst_gpu(gddIdx, :, :) = timeEst;
        gddEst_gpu(gddIdx, :, :) = gddEst;
    end
    
    fprintf('FCT computation completed: %.2f seconds\n', toc(ticFCT));
    
    %% === 4. Squeezing Compression (no data transfer) ===
    if sqNum > 0 && sqThr > 0
        fprintf('Starting GPU squeezing compression...\n');
        ticSqueeze = tic;
        
        tfComp_gpu = Msqueeze0_gpu_core(x_gpu, tf0_gpu, gddEst_gpu, timeEst_gpu, ...
                                        t, gddAxis, sqNum, sqThr);
        
        fprintf('GPU compression completed: %.2f seconds\n', toc(ticSqueeze));
    else
        tfComp_gpu = tf0_gpu;
    end
    
    %% === 5. Transfer Results Back to CPU ===
  
    ticTransfer = tic;
    
    tfComp = gather(tfComp_gpu);
    tfOrig = gather(tf0_gpu);
    
    fprintf('Transfer completed: %.2f seconds\n', toc(ticTransfer));
    
    %% === 6. Total Execution Time ===
    totalTime = toc(ticTotal);
    fprintf('Integrated GPU pipeline completed! Total time: %.2f seconds\n', totalTime);
end

function tfComp_gpu = Msqueeze0_gpu_core(x_gpu, tf0_gpu, gdd_gpu, time_gpu, tVec, gddVec, numIter, thresh)
    % MSQUEEZE0_GPU_CORE - Core squeezing algorithm on GPU
    %
    % INPUTS:
    %   x_gpu     - Input signal on GPU
    %   tf0_gpu   - Time-frequency representation on GPU
    %   gdd_gpu   - Estimated GDD matrix on GPU
    %   time_gpu  - Estimated time matrix on GPU
    %   tVec      - Time vector
    %   gddVec    - GDD vector
    %   numIter   - Number of squeezing iterations
    %   thresh    - Threshold factor
    %
    % OUTPUTS:
    %   tfComp_gpu - Compressed time-frequency representation on GPU
    
    [Ngdd, Nt, Nf] = size(tf0_gpu);  % Ngdd: GDD points, Nt: time points, Nf: frequency points
    
    % Parameters (all computations on GPU)
    Threshold = thresh * mean(abs(x_gpu).^2);
    dt = single(tVec(2) - tVec(1));
    dgdd = single(gddVec(2) - gddVec(1));
    midGdd = floor(length(gddVec)/2) + 1;
    
    % Index conversion for GDD and time (on GPU)
    gdd_idx = round(gdd_gpu ./ dgdd) + midGdd;
    time_idx = round(time_gpu ./ dt) + 1;
    gdd_idx = real(gdd_idx);
    time_idx = real(time_idx);
    
    fprintf('GPU compression core: %d iterations\n', numIter);
    
    for iter = 1:numIter
        fprintf('Iteration %d: ', iter);
        ticIter = tic;
        
        % Threshold filtering
        mask = abs(tf0_gpu) > Threshold;
        nPoints = nnz(mask);
        
        if nPoints > 0
            % Get all valid points
            [rows, cols, slices] = ind2sub([Ngdd, Nt, Nf], find(mask));
            
            % Source values
            srcVals = tf0_gpu(mask);
            
            % Target positions based on GDD and time estimates
            targetGdd = gdd_idx(mask);
            targetTime = time_idx(mask);
            
            % Validate target positions
            valid = targetGdd >= 1 & targetGdd <= Ngdd & ...
                    targetTime >= 1 & targetTime <= Nt;
            
            if any(valid)
                srcVals = srcVals(valid);
                targetGdd = targetGdd(valid);
                targetTime = targetTime(valid);
                slices = slices(valid);
                
                %% === Sparse matrix accumulation ===
                % Create linear indices
                idxDest = sub2ind([Ngdd, Nt, Nf], targetGdd, targetTime, slices);
                
                % Accumulate to new matrix
                tfNew_gpu = accumarray(idxDest, srcVals, [Ngdd*Nt*Nf, 1]);
                tfNew_gpu = reshape(tfNew_gpu, Ngdd, Nt, Nf);
            else
                tfNew_gpu = gpuArray.zeros(Ngdd, Nt, Nf, 'single');
            end
        else
            tfNew_gpu = gpuArray.zeros(Ngdd, Nt, Nf, 'single');
        end
        
        tf0_gpu = tfNew_gpu;
        fprintf('Processed %d points, %.3f seconds\n', nPoints, toc(ticIter));
    end
    
    tfComp_gpu = tf0_gpu;
end  