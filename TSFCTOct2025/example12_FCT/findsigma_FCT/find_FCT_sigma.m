
function [opt_gs1,entro1]=find_FCT_sigma(f,Hz,gs,chrrange,l)


% This function determines the optimal Gaussian window width (sigma) by minimizing
% the Renyi entropy of the time-frequency-chirp rate distribution.
%
% Inputs:
%   f         : Input signal (vector, Nx1)
%   Hz        : Sampling frequency in Hz (positive scalar)
%   gs        : Candidate sigma values vector (1xM)
%   chrrange  : Chirp rate range limit in 1/Hz² (determines GDD axis resolution)
%   l         : Order of Renyi entropy (positive scalar, typically 2 or 3)
%
% Outputs:
%   opt_gs1   : Optimal sigma value that minimizes Renyi entropy (scalar)
%   entro1    : Minimum Renyi entropy value corresponding to opt_gs1 (scalar)
%








leng=length(gs); 
entro1=zeros(1,leng); 
h = waitbar(0, 'Processing...');
for j = 1:length(gs)
   [tfc] = FCT(f,Hz,gs(j),chrrange);
    upp1 = sum(sum(sum(abs(tfc).^(2*l))));
    dow1 = sum(sum(sum(abs(tfc).^2)))^(l);
  
    
    entro1(j) = 1/(1-l)*log(upp1./dow1);
    
    % 更新进度条
    waitbar(j/length(gs), h, sprintf('Progress: %d', j, length(gs)));
    
    results = [j, gs(j), entro1(j)];
    disp(results);
end
close(h); % 关闭进度条

[gs_min1, min_ind1]=min(entro1);


opt_gs1=gs(min_ind1);
disp(opt_gs1);
