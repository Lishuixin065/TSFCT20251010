function [tfr0,Ts] = SST_2nd_hop(x,Hz,sigma,hop)
%   STFT
%	x       : Signal.
%	hlength : Window length.
% hop: hop size. The hop size equals to hlength - num overlapped points
%   tfr1  : STFT 

%ran@cau.edu.cn

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
    



tLen = length(t) ;% number of bins of time
fLen =length(tfrtic); % number of bins of frequency
ht = -tLen:tLen ;
time=L*ht/tLen/hop;
Lh=tLen;

tfr0 =zeros( fLen, tLen);
tfr1 =zeros( fLen, tLen);
tfr2 =zeros( fLen, tLen);
h=exp(-1/(2*sigma^2)*time.^2');
th=time'.*exp(-1/(2*sigma^2)*time.^2');
tth=time.^2'.*exp(-1/(2*sigma^2)*time.^2');


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
  
      omega =zeros( fLen, tLen);
      del_f=tfrtic(2)-tfrtic(1);
       omega1 =real((tfr0.*tfr1)./(tfr1.*tfr1-tfr0.*tfr2)./(2.0*pi*1i));
       omega1=round(omega1./del_f);
      
        for i=1:fLen
           omega(i,:)=i-omega1(i,:);
        end    
     
Ts=zeros(fLen,tLen);
    for b=1:tLen%time
    % Reassignment step
    for eta=1:fLen%frequency
        if abs(tfr0(eta,b))>0.0001%you can set much lower value than this.
            k = omega(eta,b)+1;
            if k>=1 && k<=fLen
                Ts(k,b) = Ts(k,b) + tfr0(eta,b);
            end
        end
    end
    end  
  Ts=Ts/(tLen);    
   