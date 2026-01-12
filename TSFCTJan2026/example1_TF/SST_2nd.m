function [Ts,tfr, tfrtic,omega] = SST_2nd(x, Hz, s)


	%% check input signals
 [xrow,~] = size(x);
  if (xrow~=1)
        x = x';
  end
  N = length(x);
   if mod(N,2)==0
        f = (0:N/2)*Hz/N;
    else
        f = (0:(N-1)/2)*Hz/N;
   end
   tfrtic =f ;


Lh=N;
ht = -N:N ;
time=ht/Hz;
tLen = length(x) ;% number of bins of time
fLen =length(tfrtic) ; % number of bins of frequency

   g0 = @(t) exp(-2*pi^2*s^2*t.^2);
   tg0 = @(t) exp(-2*pi^2*s^2*t.^2).*t; 
   ttg0 = @(t) exp(-2*pi^2*s^2*t.^2).*t.^2;
   

h0=1;
	%% run STFT and reassignment rule
tfr =zeros( fLen, tLen); 
tfr0 =zeros( fLen, tLen); 
tfr1 =zeros( fLen, tLen); 
tfr2 =zeros( fLen, tLen); 
omega=zeros( fLen, tLen);
Ts=zeros( fLen, tLen);
    for tidx = 1:tLen      
        ti = tidx;
        tau = -min([round(N/2)-1,ti-1]):min([round(N/2)-1,N-ti]);
        indices= rem(N+tau,N)+1;      
        tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ; tf2 = zeros(N, 1) ;
        g=g0(time(Lh+1+tau));
        tg=tg0(time(Lh+1+tau));
        ttg=ttg0(time(Lh+1+tau));
        tf0(indices) = x(ti+tau).*g;
        tf1(indices) = x(ti+tau).*tg;
        tf2(indices) = x(ti+tau).*ttg;
        tf0 = fft(tf0) ; tf0 = tf0(1:fLen) ;
        tf1 = fft(tf1) ; tf1 = tf1(1:fLen) ;
        tf2 = fft(tf2) ; tf2 = tf2(1:fLen) ;
        tfr(:, tidx) = tf0(1:fLen) ;   
        tfr0(:, tidx) = tf0(1:fLen);  
        tfr1(:, tidx) = tf1(1:fLen);
        tfr2(:, tidx) = tf2(1:fLen);
             
    end
       del_f=tfrtic(2)-tfrtic(1); 
       omega1 =real((tfr0.*tfr1)./(tfr1.*tfr1-tfr0.*tfr2)./(2.0*pi*1i));
       omega1=round(omega1./del_f);

        for i=1:fLen
           omega(i,:)=i-omega1(i,:);
        end    


    for b=1:tLen%time
    for eta=1:fLen%frequency
        if abs(tfr(eta,b))>eps%you can set much lower value than this.
            k = omega(eta,b);
            if k>=1 && k<=fLen
                Ts(k,b) = Ts(k,b) + tfr(eta,b);
            end
        end
    end
    end  

    Ts=Ts/N/h0;
end







