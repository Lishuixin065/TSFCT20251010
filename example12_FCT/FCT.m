function [tfc0] = FCT(x,fs,s,chrrange)
% % Inputs:
%   x         : Input signal (vector)
%   fs        : Sampling frequency (Hz)
%   s         : Standard deviation of Gaussian window (in samples)
%   chrrange  : Maximum chirp rate range (determines GDD axis resolution)
%
% Outputs:
%   tfc0      : FCT coefficients (using base Gaussian window g)


    %% 判断是否为列向量，是则转置
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    %% 预处理
    N = length(x);
    tao = (0:N-1)/fs;   
    t = (0:N-1)/fs; 
    if mod(N,2)==0
        f = (0:N/2)*fs/N;
    else
        f = (0:(N-1)/2)*fs/N;
    end
    L = length(f);
  

   crate=linspace(-chrrange,chrrange,2*round(N/2)+1);
   tfrtic =f ;
   tcrtic =crate;
   tLen = length(t(1:length(x))) ;
   fLen =length(tfrtic) ;
   cLen = length(crate); 
  
   	%% run STFT and reassignment rule
   tfc0 =zeros(cLen, tLen,fLen);
 
    
    %% 计算
 fprintf(['ahirp-rate total: ',num2str(cLen), '; now:     ']) ;
    for cidx = 1:cLen
    fprintf('\b\b\b\b') ;	
    tmp = sprintf('%4d',cidx) ; 
    fprintf(tmp) ;
    ahirp = crate(cidx);
   

    gt0 =FCT_windows(s,ahirp,0);
       for tidx = 1:N
        gh0 = gt0(t(tidx)-tao); 
        gh0= conj(gh0);
        tf0 = fft(gh0 .* x);
        tf0=tf0(1:L);
        tfc0(cidx, tidx,:) = tf0(1:L); 
        end          

    end
fprintf('\n') ; 


    

         
end
