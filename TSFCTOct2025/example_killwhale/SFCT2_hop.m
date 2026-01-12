function [tfc0,tfrtic,tcrtic,T,R] = SFCT2_hop(x,fs,s,chrrange,hop)
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
      t = (0:hop:N-1)/fs; 
     if mod(N,2)==0
        f_full = (0:round(N/2))*fs/N;
    else
        f_full = (0:round((N-1)/2))*fs/N;
    end
  
   f = f_full(1:hop:end);             % 最后才下采样  

   
   
    tfrtic =f ;
    tLen = length(t);
   fLen =length(tfrtic) ;
   crate=linspace(-chrrange,chrrange,2*round(tLen/4)+1); 
   tcrtic =crate;
   cLen = length(crate); 
 
   	%% run STFT and reassignment rule
   tfc0 =zeros(cLen, tLen,fLen);
   tfc1 =zeros(cLen, tLen,fLen);
   tfc2 =zeros(cLen, tLen,fLen);
   T =zeros(cLen, tLen,fLen);
   R =zeros(cLen, tLen,fLen);

   del_t=hop/fs;
   

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
    for tidx = 1:tLen
        gh0 = gt0(t(tidx)-tao); 
        gh0 = conj(gh0);
        tf0_full = fft(gh0 .* x);
        tf0_full = tf0_full(1:length(f_full));  
        tf0 = tf0_full(1:hop:end); 
      
       

        gh1 = gt1(t(tidx)-tao);
        gh1 = conj(gh1);
        tf1_full = fft(gh1 .* x);
        tf1_full = tf1_full(1:length(f_full));  
        tf1 = tf1_full(1:hop:end);            
        tfc1(cidx, tidx,:) = tf1;
       

        gh2 = gt2(t(tidx)-tao);
        gh2 = conj(gh2);     
        tf2_full = fft(gh2 .* x);
        tf2_full = tf2_full(1:length(f_full));  
        tf2 = tf2_full(1:hop:end);            
        tfc2(cidx, tidx,:) = tf2;
       

              

     M0=tf2.*tf0-tf1.*tf1;
     M1=-tf0.*tf1;
     M2 =(tf0.*tf0);    
    
       
             
    
     gamma=eps;
     tf0(abs(tf0)<gamma & abs(M0) < gamma)=0; 
     tfc0(cidx, tidx,:) = tf0;
     omega1=(tidx-1)*del_t+M1./M0./(2.0*pi*1i);
     lambda1 =(ahirp)+M2./M0./(2.0*pi*1i); 

     omega1(abs(tf0)<gamma & abs(M0) < gamma ) =Inf;
     lambda1(abs(tf0)<gamma & abs(M0) < gamma ) =Inf;
     
        


     R(cidx, tidx,:)=lambda1;
     T(cidx,tidx,:)=omega1;  


      
     
        end          

    end
fprintf('\n') ; 

         
end
