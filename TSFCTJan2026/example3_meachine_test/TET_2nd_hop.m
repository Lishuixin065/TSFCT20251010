function[Ts,tfr0,t,f, GroupDelay] = TET_2nd_hop_hop(x,Hz,sigma,hop)

%   STFT
%	x       : Signal.
%	hlength : Window length.
% hop: hop size. The hop size equals to hlength - num overlapped points
%   tfr1  : STFT 



[xrow,xcol] = size(x) ;
	%% check input signals
if (xcol~=1)
    error('x must have only one column');

end


[xrow,xcol] = size(x);
N=xrow;
t=1:hop:N;
L=N/Hz;
%tfrtic = (0:hop:round(N/2))/L;
  if mod(N,2)==0
        tfrtic = (0:hop:round(N/2))*Hz/N;
  else
        tfrtic = (0:hop:round((N-1)/2))*Hz/N;
  end
    
f=tfrtic;


tLen = length(t) ;% number of bins of time
fLen =length(tfrtic); % number of bins of frequency
ht = -tLen:tLen ;
time=L*ht/tLen/hop;
Lh=tLen;

tfr0 =zeros( fLen, tLen);
tfr1 =zeros( fLen, tLen);
tfr2 =zeros( fLen, tLen);
h=exp(-(2*pi^2*sigma^2)*time.^2');
th=time'.*exp(-(2*pi^2*sigma^2)*time.^2');
tth=time.^2'.*exp(-(2*pi^2*sigma^2)*time.^2');


 for tidx = 1:tLen      
        ti = t(tidx);
        tau=-min([round(tLen/2)-1,Lh,ti-1]):min([round(tLen/2)-1,Lh,xrow-ti]);
        indices= rem(tLen+tau,tLen)+1;      
        tf0 = zeros(tLen, 1); tf1 = zeros(tLen, 1);  tf2 = zeros(tLen, 1);  
        tf0(indices)=x(ti+tau).*h(Lh+1+tau);
        tf1(indices)=x(ti+tau).*th(Lh+1+tau);
        tf2(indices)=x(ti+tau).*tth(Lh+1+tau);
        tf0 = fft(tf0); tf1 = fft(tf1); tf2 = fft(tf2); 
        tf0 = tf0(1:fLen); tf1 = tf1(1:fLen); tf2 = tf2(1:fLen);
        tfr0(:, tidx) = tf0(1:fLen) ; 
        tfr1(:, tidx) = tf1(1:fLen) ; 
        tfr2(:, tidx) = tf2(1:fLen) ; 
  end
  
t0=t/Hz;del_t=hop/Hz;
      for fidx=1:length(f)
        for tidx=1:length(t)
          tfr0(fidx,tidx)=tfr0(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t0(tidx));
          tfr1(fidx,tidx)=tfr1(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t0(tidx));
          tfr2(fidx,tidx)=tfr2(fidx,tidx).*exp(-2*pi*1i*f(fidx)*t0(tidx));
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
           if m==b
            Ts(prt,m) = Ts(prt,m) + tfr0(prt,b);
           end
        end
    end

    

         
end
