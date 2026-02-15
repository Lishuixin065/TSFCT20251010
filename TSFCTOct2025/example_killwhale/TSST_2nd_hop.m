function[Tx,tfc0,t,f,T] = TSST_2nd_hop(x,fs,s,hop)

   [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    %% 预处理
    N = length(x);

  
    tao = (0:N-1)/fs;   
    t = (0:hop:N-1)/fs; 
     if mod(N,2)==0
        f_full = (0:round(N/2))*fs/N;
    else
        f_full = (0:round((N-1)/2))*fs/N;
    end
  
   f = f_full(1:hop:end);             % 最后才下采样 
   del_t=t(2)-t(1);
   tfrtic =f ;
   tLen = length(t);
   fLen =length(tfrtic) ;
   	%% run STFT and reassignment rule
   tfc0 =zeros( tLen,fLen);
   tfc1 =zeros( tLen,fLen);
   tfc2 =zeros(tLen,fLen);
   T =zeros( tLen,fLen);
    

 
    gt0 =FCT_windows(s,0,0);
    gt1 = FCT_windows(s,0,1);
    gt2 = FCT_windows(s,0,2);
    for tidx = 1:tLen
        gh0 = gt0(t(tidx)-tao); 
        gh0 = conj(gh0);  
        tf0_full = fft(gh0 .* x);
        tf0_full = tf0_full(1:length(f_full));  % 正频率
        tf0 = tf0_full(1:hop:end);              % 频率下采样   
        tfc0(tidx, :) = tf0;

        gh1 = gt1(t(tidx)-tao);
        gh1 = conj(gh1);
        tf1_full = fft(gh1 .* x);
        tf1_full=tf1_full(1:length(f_full));
        tf1 = tf1_full(1:hop:end);
        tfc1(tidx,:) = tf1; 

         gh2 = gt2(t(tidx)-tao);
         gh2 = conj(gh2);
         tf2_full = fft(gh2 .* x);
         tf2_full=tf2_full(1:length(f_full));
         tf2 = tf2_full(1:hop:end);
         tfc2(tidx,:) = tf2; 
         


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
           if 1<=m && m<=tLen
            Tx( m,prt) = Tx( m,prt) + tfc0(b,prt);
           end
        end
    end


    

         
end
