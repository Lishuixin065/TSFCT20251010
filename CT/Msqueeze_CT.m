
%% Algorithm for multiple compressions
function [tfrsq2] = Msqueeze_CT(x, tfrsq1, lambda, omega, tfrtic,tcrtic, num)


[a, b, c] = size(tfrsq1);

Threshold=0.0001*mean(abs(x).^2);
 del_c=tcrtic(2)-tcrtic(1);
 del_f=tfrtic(2)-tfrtic(1);

omega=round(omega./del_f)+1;

for cidx=1:a
    for tidx=1:c     
 
      lambda(cidx,:,tidx)=round((lambda(cidx,:,tidx)+tcrtic(end))./del_c)+1;
    
    end
end



valid_indices = (tfrsq1 < Threshold) & (omega < Threshold);
omega(valid_indices)=0; 
lambda(valid_indices)=0;










fprintf(['Total time iterations: ', num2str(c), '; Current:     ']);

for kk = 1 : num
    tfrsq2 = zeros(a, b, c);
    for tidx = 1 : c
        fprintf('\b\b\b\b');  % Backspace to update the progress
        tmp = sprintf('%4d', tidx);
        fprintf(tmp);  % Print the current iteration index
        lambda2 = squeeze(lambda(:, :, tidx));
        omega2 = squeeze(omega(:,:,tidx));
        tf0 = squeeze(tfrsq1(:,:,tidx));
        sst = zeros(a, b);
      
        
        for cidx = 1 : a  % Loop over chirps
            for fidx = 1 : b  % Loop over frequencies
                k = omega2(cidx, fidx); m = lambda2(cidx, fidx);
                if abs(tf0(cidx, fidx)) > Threshold  % Filter based on the threshold
                    if (k <=b ) && (k >=1) ...
                       && (m >= 1) && (m <= a)
                        sst(m , k ) = sst(m, k ) + tf0(cidx, fidx);
                    end
                end
            end
        end
        tfrsq2(:,:,tidx) = sst;
    end
    tfrsq1 = tfrsq2;  % Update the initial TFR with the compressed data
end

fprintf('\n');  % New line after the progress is completed
end