function [Ts,tfr, tfrtic] = TSST_2nd(x,Hz,sigma)

    %% 判断是否为列向量，是则转置
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
    %% 预处理
    N = length(x);
    E = mean(abs(x)); 
    if mod(N,2)==0
        f = (0:N/2)*Hz/N;
    else
        f = (0:(N-1)/2)*Hz/N;
    end
    L = length(f);
    del_t = 1/Hz;
    tfrtic=f;



   Lh=N;
   ht = -N:N ;
   time=ht/Hz;
   tLen = length(x) ;% number of bins of time
   fLen =length(tfrtic) ; % number of bins of frequency







    %% 计算
    %计算Wx 
    gt0 = @(t) sigma^(-1)*(2*pi)^(-1/2).*exp(-t.^2/sigma^2/2);
    gt1 = @(t) sigma^(-1)*(2*pi)^(-1/2).*exp(-t.^2/sigma^2/2).*t; 
    gt2 = @(t) sigma^(-1)*(2*pi)^(-1/2).*exp(-t.^2/sigma^2/2).*t.^2;
   
     tfr = zeros(fLen, tLen); 
    tfr0 = zeros(fLen, tLen);   
    tfr1 = zeros(fLen, tLen);
    tfr2 = zeros(fLen, tLen); 


    for tidx = 1:tLen      
        ti = tidx;
        tau = -min([round(N/2)-1,ti-1]):min([round(N/2)-1,N-ti]);
        indices= rem(N+tau,N)+1;      
        tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ; tf2 = zeros(N, 1) ;
        g=gt0(time(Lh+1+tau));
        tg=gt1(time(Lh+1+tau));
        ttg=gt2(time(Lh+1+tau));
        tf0(indices) = x(ti+tau).*g;
        tf1(indices) = x(ti+tau).*tg;
        tf2(indices) = x(ti+tau).*ttg;
        tf0 = fft(tf0) ; tf0 = tf0(1:fLen) ;
        tf1 = fft(tf1) ; tf1 = tf1(1:fLen) ;
        tf2 = fft(tf2) ; tf2 = tf2(1:fLen) ;
        tfr(:, tidx) = tf0(1:fLen) ;   
        tfr0(:, tidx) = tf0(1:fLen) ;  
        tfr1(:, tidx) = tf1(1:fLen) ;
        tfr2(:, tidx) = tf2(1:fLen) ;     
    end

    t = (0:N-1)/Hz; 
      for fidx=1:length(f)
        for tidx=1:length(t)
          tfr0(fidx,tidx)=tfr0(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t(tidx));
          tfr1(fidx,tidx)=tfr1(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t(tidx));
          tfr2(fidx,tidx)=tfr2(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t(tidx));
        end
       end



tfr1=-tfr1./(2*pi*1i)/sigma^2;
tfr2=(-tfr0/sigma^2+tfr2/sigma^4)*(-1/(2*pi*1i))^2;

        %% 求群延迟
            Denominator = tfr2.*tfr0-tfr1.*tfr1;%似乎严格按照（25）式来的
            Numerator = tfr0.*tfr1;
            Numerator2 =tfr0.*tfr0;
            p = Numerator./Denominator./(2*pi*1i);
            q = Numerator2./Denominator./(2*pi*1i);
            for tau = 1:tLen
                p(:,tau) =(tau)*del_t+p(:,tau);
            end
             
          GroupDelay =round(real(p)./del_t);%群延迟估计 

            

        
        %% TFR
     Ts = zeros(fLen,tLen);
    for prt=1:fLen
        for b=1:tLen
           m= GroupDelay(prt,b);
           if 1<=m&& m<=N
            Ts(prt,m) = Ts(prt,m) + tfr0(prt,b);
           end
        end
    end

end




