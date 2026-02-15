function [tfc0,tfc1,tfc2,tfrtic,tcrtic,T,R] = SFCT2(x,fs,s,chrrange)
% Inputs:
%   x         : Input signal (vector)
%   fs        : Sampling frequency (Hz)
%   s         : Standard deviation of Gaussian window (in samples)
%   chrrange  : Maximum chirp rate range (determines GDD axis resolution)
%
% Outputs:
%   tfc0      : FCT coefficients (using base Gaussian window g)
%   tfc1      : FCT coefficients (using first-order window tg)
%   tfc2      : FCT coefficients (using second-order window t²g)
%   tfrtic    : Frequency axis (Hz)
%   tcrtic    : Group Delay Dispersion (GDD) axis 
%   T         : Second order GD reassignment operator 
%   R         : Second order  GDD reassignment operator 



    %% 判断是否为列向量，是则转置
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    %% 预处理
    N = length(x);
    E = mean(abs(x));
    tao = (0:N-1)/fs;   
    t = (0:N-1)/fs; 
    if mod(N,2)==0
        f = (0:round(N/2))*fs/N;
    else
        f = (0:round((N-1)/2))*fs/N;
    end
    

    L = length(f);
    dt = 1/fs;

   crate=linspace(-chrrange,chrrange,2*round(N/2)+1);
   tfrtic =f ;
   tcrtic =crate;
   tLen = length(t(1:length(x))) ;
   fLen =length(tfrtic) ;
   cLen = length(crate); 
  del_t=t(2)-t(1);
  del_c=tcrtic(2)-tcrtic(1);
   	%% run STFT and reassignment rule
   tfc0 =zeros(cLen, tLen,fLen);
   tfc1 =zeros(cLen, tLen,fLen);
   tfc2 =zeros(cLen, tLen,fLen);
   T =zeros(cLen, tLen,fLen);
   R =zeros(cLen, tLen,fLen);
          del_t=1/fs;
     del_c=tcrtic(2)-tcrtic(1);
     del_f= tfrtic(2)-tfrtic(1);

    %% 计算
 fprintf(['ahirp-rate total: ',num2str(cLen), '; now:     ']) ;
    for cidx = 1:cLen
    fprintf('\b\b\b\b') ;	
    tmp = sprintf('%4d',cidx) ; 
    fprintf(tmp) ;
    ahirp = crate(cidx);
   

    gt0 =FCT_windows(s,ahirp,0);
    gt1 = FCT_windows(s,ahirp,1);
    gt2 = FCT_windows(s,ahirp,2);
    for tidx = 1:N
        gh0 = gt0(t(tidx)-tao); 
        gh0 = conj(gh0);
        tf0 = fft(gh0 .* x);
        tf0=tf0(1:L);
        tfc0(cidx, tidx,:) = tf0(1:L);
 
        gh1 = gt1(t(tidx)-tao);
        gh1 = conj(gh1);
        tf1 = fft(gh1 .* x);
        tf1=tf1(1:L);
        tfc1(cidx, tidx,:) = tf1(1:L); 

         gh2 = gt2(t(tidx)-tao);
         gh2 = conj(gh2);
         tf2 = fft(gh2 .* x);
         tf2=tf2(1:L);
         tfc2(cidx, tidx,:) = tf2(1:L);                

     M0=tf2.*tf0-tf1.*tf1;
     M1=-tf0.*tf1;
     M2 =(tf0.*tf0);    
    
       
     omega1=(tidx-1)*del_t+M1./M0./(2.0*pi*1i);
     lambda1 =(ahirp)+M2./M0./(2.0*pi*1i); 
     R(cidx, tidx,:)=lambda1;
     T(cidx,tidx,:)=omega1;  


      
     
        end          

    end
fprintf('\n') ; 


    

         
end
