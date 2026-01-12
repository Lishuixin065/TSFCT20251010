function[tfc0,tfrtic] = STFT(x,fs,s)

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
        f = (0:N/2)*fs/N;
    else
        f = (0:(N-1)/2)*fs/N;
    end
  
   
   L = length(f);
   tfrtic =f ;
   tLen = length(t(1:length(x))) ;
   fLen =length(tfrtic) ;
   	%% run STFT and reassignment rule
   tfc0 =zeros( tLen,fLen);

    


     gt0=@(t) exp(-2*pi^2*s^2*t.^2);
 

       for tidx = 1:N
        gh0 = gt0(t(tidx)-tao); 
        gh0 = conj(gh0);
        tf0 = fft(gh0 .* x);
        tf0=tf0(1:L);
        tfc0( tidx,:) = tf0(1:L);   
        end          
         
end
