function [tfc0,tfc1,tfc2,tfrtic,aa,R,T] = FCT2nd_hop(x,Hz,sigma,arange,hop)

%   STFT
%	x       : Signal.
%	hlength : Window length.
% hop: hop size. The hop size equals to hlength - num overlapped points
%   tfr1  : STFT 



[~,xcol] = size(x) ;
	%% check input signals
if (xcol~=1)
    error('x must have only one column');

end

N=length(x);
t0=0:hop:N-1;
t=t0/Hz;
del_t=hop/Hz;

%tfrtic = (0:hop:round(N/2))/L;
  if mod(N,2)==0
        tfrtic = (0:hop:N/2)*Hz/N;
  else
        tfrtic = (0:hop:(N-1)/2)*Hz/N;
  end

f=tfrtic;


tLen = length(t) ;% number of bins of time
fLen =length(tfrtic); % number of bins of frequency
ht = -tLen:tLen ;
Lh=tLen;
time=ht/Hz;

aa= linspace(-arange, arange, 2*round(tLen/2)+1);  % ahirp rate vector
aLen=length(aa)  ;



tfc0 =zeros(aLen,tLen,fLen);
tfc1 =zeros(aLen,tLen,fLen);
tfc2 =zeros(aLen,tLen,fLen);
 T=zeros(aLen,tLen,fLen);
 R=zeros(aLen,tLen,fLen);


 for cidx = 1:aLen
        fprintf('\b\b\b\b');
        fprintf('%4d', cidx);
        
        ahirp = aa(cidx);  % current ahirp rate
    gt0 =FCT_windows(sigma,ahirp,0);
    gt1 = FCT_windows(sigma,ahirp,1);
    gt2 = FCT_windows(sigma,ahirp,2);

    h = conj(gt0(time'));
    th=conj(gt1(time'));
    tth=conj(gt2(time'));




 for tidx = 1:tLen      
        ti = t0(tidx)+1;
        tau=-min([round(tLen/2)-1,tLen,ti-1]):min([round(tLen/2)-1,tLen,N-ti]);    
        indices= rem(tLen+tau,tLen)+1;      
        tf0 = zeros(tLen, 1); tf1 = zeros(tLen, 1);  tf2 = zeros(tLen, 1);  
        %tf0(indices)=x(ti+tau).*h(Lh+1+tau);
        tf0(indices)=x(ti+tau).*h(Lh+1-tau);
        tf1(indices)=x(ti+tau).*th(Lh+1-tau);
        tf2(indices)=x(ti+tau).*tth(Lh+1-tau);

        % tf0(indices)=x(ti+tau).*conj(gt0(time(tLen+1-tau)))';  
        % tf1(indices)=x(ti+tau).*conj(gt1(time(tLen+1-tau)))';
        % tf2(indices)=x(ti+tau).*conj(gt2(time(tLen+1-tau)))';




        tf0 = fft(tf0); tf1 = fft(tf1); tf2 = fft(tf2); 
        tf0 = tf0(1:fLen); 
        tf1 = tf1(1:fLen);
        tf2 = tf2(1:fLen);
         
        
        tf0=tf0.*exp(-2*pi*1i*f'*t(tidx));
        tf1=tf1.*exp(-2*pi*1i*f'*t(tidx));
        tf2=tf2.*exp(-2*pi*1i*f'*t(tidx));


        tfc0(cidx, tidx,:) = tf0; 
        tfc1(cidx, tidx,:) = tf1 ; 
        tfc2(cidx,tidx,:) = tf2; 

    M0=tf2.*tf0-tf1.*tf1;
     M1=-tf0.*tf1;
     M2 =(tf0.*tf0);    
    
     thre=0.001*mean(abs(M0));   
     E0 = abs(M0) > thre; 

     omega1=(tidx-1)*del_t+M1./M0./(2.0*pi*1i);
     lambda1 =(ahirp)+M2./M0./(2.0*pi*1i); 
     R(cidx, tidx,:)=lambda1.*E0;
     T(cidx,tidx,:)=omega1.*E0;  


  end
  
 end
 fprintf('\n');
end