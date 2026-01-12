function [ct, lam, omega, f, c] = CT_2nd(x, fs, sigma, crange)

    % Input renamed:
    % x -> signal
    % Hz -> fs (sampling frequency)
    % sigma -> sigma (Gaussian window parameter)
    % chrrange -> crange (chirp rate range)
    
    % Ensure row vector
    [nrow, ~] = size(x);
    if nrow ~= 1
        x = x';
    end

    % Signal parameters
    N = length(x);
    
    % Chirp rate vector
    c = linspace(-crange, crange, 2*round(N/2)+1);  % chirp rate vector
    
    % Window time indices and vector
    window_indices = -N:N;  % indices for window
    window_time = window_indices / fs;  % time vector for window
    
    L = 1/sigma/sqrt(2*pi);  % normalization factor
    
    % Window function handles
    g = @(time, lam) L .* exp(-0.5*time.^2/sigma^2) .* exp(-pi*1i*lam*time.^2);
    tg = @(time, lam) L .* exp(-0.5*time.^2/sigma^2) .* exp(-pi*1i*lam*time.^2) .* time;
    ttg = @(time, lam) L .* exp(-0.5*time.^2/sigma^2) .* exp(-pi*1i*lam*time.^2) .* time.^2;

    % Frequency vector
    if mod(N, 2) == 0
        f = (0:N/2) * fs/N;
    else
        f = (0:(N-1)/2) * fs/N;
    end

    % Dimensions
    ntime = N;          % time bins
    nfreq = length(f);  % frequency bins
    nchirp = length(c); % chirp rate bins

    % Initialize results
    ct = zeros(nchirp, nfreq, ntime);  % chirplet transform
    lam = zeros(nchirp, nfreq, ntime); % estimated chirp rate
    omega = zeros(nchirp, nfreq, ntime); % estimated frequency

    fprintf('Chirp-rate total: %d; now:     ', nchirp);

    for cidx = 1:nchirp
        fprintf('\b\b\b\b');
        fprintf('%4d', cidx);
        
        chirp = c(cidx);  % current chirp rate
        
        for tidx = 1:ntime
            % Window indices
            center_idx = tidx;
            tau = -min([round(N/2)-1, N, center_idx-1]):min([round(N/2)-1, N, N-center_idx]);
            
            % Circular indices
            idx = mod(N + tau, N) + 1;
            
            % Window values at current time offset
            g_vals = g(window_time(N+1+tau), chirp);
            tg_vals = tg(window_time(N+1+tau), chirp);
            ttg_vals = ttg(window_time(N+1+tau), chirp);

            % Initialize transforms
            tf0 = zeros(N, 1);
            tf1 = zeros(N, 1);
            tf2 = zeros(N, 1);
            
            % Compute chirplet transforms
            tf0(idx) = x(center_idx+tau) .* g_vals;
            tf1(idx) = x(center_idx+tau) .* tg_vals;
            tf2(idx) = x(center_idx+tau) .* ttg_vals;
            
            % FFT
            tf0 = fft(tf0)/fs;
            tf0 = tf0(1:nfreq);
            tf1 = fft(tf1)/fs;
            tf1 = tf1(1:nfreq);
            tf2 = fft(tf2)/fs;
            tf2 = tf2(1:nfreq);
           
            % Step sizes
            dc = c(2) - c(1);  % chirp rate step
            df = f(2) - f(1);  % frequency step

            % Compute moments for reassignment
            M0 = tf0 .* tf2 - tf1 .* tf1;
            M1 = tf0 .* tf1;
            M2 = -(tf0 .* tf0);
            
            % Threshold small values
            tf0(abs(M0) < 1e-10) = 0;
            ct(cidx, :, tidx) = tf0;

            % Reassignment
            omega1 = imag(M1./M0/(2*pi));
            lam1 = imag(M2./M0/(2*pi));
            
            lam2 = chirp + lam1;
            omega2 = zeros(size(omega1));
            
            % Frequency correction
            for i = 1:length(omega1)
                omega2(i) = (i-1)*df + omega1(i);
            end
             
            % Store estimates
            lam(cidx, :, tidx) = lam2;
            omega(cidx, :, tidx) = omega2;
        end
    end
    
    fprintf('\n');
end