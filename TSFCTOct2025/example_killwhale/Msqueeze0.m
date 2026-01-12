
%% Algorithm for multiple compressions
function [tfrsq1] = Msqueeze0(x, tfc0, R0, T0,t,tcrtic, num,thre)
%
% Inputs:
%   x       : Input signal (vector, Nx1)
%   tfc0    : Initial Frequency-domian Chirplet Transform (FCT) 
%   R0      : Group Delay Dispersion (GDD) reassignment operator
%   T0      : Group Delay (GD) reassignment operator 
%   t       : Time axis vector (1xN)
%   tcrtic  : GDD axis vector 
%   num     : Number of squeezing iterations (positive integer)
%   thre    : Energy concentration threshold (0 < thre < 1)
%
% Output:
%   tfrsq1  : HTSFCT
%

[a, b, c] = size(tfc0); 
Threshold=thre* mean(abs(x).^2);
  R=zeros(a,b,c);
  T=zeros(a,b,c);
     del_t=t(2)-t(1);
     del_c=tcrtic(2)-tcrtic(1);
     R0 =real(R0);
     T0 =real(T0);
for cidx=1:a
    for tidx=1:b 
     lam1=round(R0(cidx,tidx,:)./del_c)+floor(length(tcrtic)/2)+1; 
     R(cidx,tidx,:)=lam1;      
     t1=round(T0(cidx,tidx,:)./del_t)+1;   
     T(cidx,tidx,:)=t1;  
    end
end




fprintf(['Total time iterations: ', num2str(c), '; Current:     ']);

for kk = 1 : num
    tfrsq1 = zeros(a, b, c);
    for fidx = 1 : c
        fprintf('\b\b\b\b');  % Backspace to update the progress
        tmp = sprintf('%4d', fidx);
        fprintf(tmp);  % Print the current iteration index
        R2 = squeeze(R(:, :, fidx));
        T2 = squeeze(T(:,:,fidx));
        tf0 = squeeze(tfc0(:,:,fidx));
        sst = zeros(a, b);
   
        
        for cidx = 1 : a  % Loop over chirps
            for tidx = 1 : b  % Loop over frequencies
                k = T2(cidx, tidx); m = R2(cidx, tidx);
                if abs(tf0(cidx, tidx)) > Threshold  % Filter based on the threshold
                    if (k <=b) && (k >=1)&& (m >= 1) && (m <= a)
                        sst(m , k ) = sst(m, k) + tf0(cidx, tidx);
                    end
                end
            end
        end
        tfrsq1(:,:,fidx) = sst;
    end
    tfc0 = tfrsq1;  % Update the initial TFR with the compressed data
end

fprintf('\n');  % New line after the progress is completed
end