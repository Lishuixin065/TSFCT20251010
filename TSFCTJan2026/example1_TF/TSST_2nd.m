function[Tx,tfc0,t,f,T] = TSST_2nd(x,fs,s)

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
   del_t=t(2)-t(1);
   tfrtic =f ;
   tLen = length(t(1:length(x))) ;
   fLen =length(tfrtic) ;
   	%% run STFT and reassignment rule
   tfc0 =zeros( tLen,fLen);
   tfc1 =zeros( tLen,fLen);
   tfc2 =zeros(tLen,fLen);
   T =zeros( tLen,fLen);
    

 
    gt0 =FCT_windows(s,0,0);
    gt1 = FCT_windows(s,0,1);
    gt2 = FCT_windows(s,0,2);
    for tidx = 1:N
        gh0 = gt0(t(tidx)-tao); 
        gh0 = conj(gh0);
        tf0 = fft(gh0 .* x);
        tf0=tf0(1:L);
        tfc0( tidx,:) = tf0(1:L);
 
        gh1 = gt1(t(tidx)-tao);
        gh1 = conj(gh1);
        tf1 = fft(gh1 .* x);
        tf1=tf1(1:L);
        tfc1(tidx,:) = tf1(1:L); 

         gh2 = gt2(t(tidx)-tao);
         gh2 = conj(gh2);
         tf2 = fft(gh2 .* x);
         tf2=tf2(1:L);
         tfc2( tidx,:) = tf2(1:L);                



     M0=tf2.*tf0-tf1.*tf1;
     M1=-tf0.*tf1;    
     Groupdelay=(tidx-1)*del_t+ M1./M0./(2.0*pi*1i);
    
   

      gamma=eps;
     tfc0(abs(tf0)<gamma & abs(M0) < gamma)=0;          
     Groupdelay(abs(tf0)<gamma & abs(M0) < gamma ) = 0;
     T(tidx,:)=round(real(Groupdelay)./del_t); 
        end          


    Tx = zeros(tLen,fLen);
    for prt=1:fLen
        for b=1:tLen
           m= T(b,prt);
           if 1<=m && m<=N
            Tx( m,prt) = Tx( m,prt) + tfc0(b,prt);
           end
        end
    end


    

         
end
