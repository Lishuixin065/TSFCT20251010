function [tfc0,lambda,omega, tfrtic, tcrtic] = CT_2nd(x, Hz, sigma,chrrange)




[xrow,xcol] = size(x) ;
	%% check input signals
if (xcol~=1)
    error('x must have only one column');

end
Lh=xrow;N=xrow;
crate=linspace(-chrrange,chrrange,2*round(N/2)+1);
tcrtic = crate;% discretization of chirp rate
%crate=linspace(-chrrange,chrrange,500);

% for tfrsq

ht = -Lh:Lh ;
time=ht/Hz;
L=1/sigma/sqrt(2*pi);
g=L.*(exp(-0.5*time.^2'/sigma^2));

t = 1:length(x) ;
 if mod(N,2)==0
        tfrtic = (0:N/2)*Hz/N;
    else
        tfrtic = (0:(N-1)/2)*Hz/N;
  end



tLen = length(t(1:length(x))) ;% number of bins of time
fLen =length(tfrtic) ; % number of bins of frequency
cLen = length(crate); % number of bins of chirp rate


	%% run STFT and reassignment rule
tfc0=zeros(cLen, fLen, tLen); % synchrosqueezed chirplet transform

lambda=zeros(cLen, fLen, tLen);omega=zeros(cLen, fLen, tLen);

fprintf(['Chirp-rate total: ',num2str(cLen), '; now:     ']) ;

for cidx = 1:cLen
    fprintf('\b\b\b\b') ;	
    tmp = sprintf('%4d',cidx) ; 
    fprintf(tmp) ;
    chirp = crate(cidx);
    for tidx = 1:tLen             
        ti = t(tidx);     
        tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);   
        indices= rem(N+tau,N)+1;      
        tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ; tf2 = zeros(N, 1) ;
        tf0(indices) = x(ti+tau).*g(Lh+1+tau).*exp(-pi*1i*chirp.*(time(Lh+1+tau)').^2); % for CT with window g
        tf1(indices) = x(ti+tau).* time(Lh+1+tau)'.* g(Lh+1+tau).*exp(-pi*1i*chirp.*(time(Lh+1+tau)').^2); % for CT with window g
        tf2(indices) = x(ti+tau).*(time(Lh+1+tau)').^2.*g(Lh+1+tau).*exp(-pi*1i*chirp.*(time(Lh+1+tau)').^2); % for CT with window g
        
        
        tf0 = fft(tf0)/Hz ; tf0 = tf0(1:fLen) ;
        tf1 = fft(tf1)/Hz ; tf1 = tf1(1:fLen) ;
        tf2 = fft(tf2)/Hz ; tf2 = tf2(1:fLen) ;
       
     
       
      
        del_c=tcrtic(2)-tcrtic(1);
        del_f=tfrtic(2)-tfrtic(1);


        M0=tf0.*tf2-tf1.*tf1   ;
        M1=(tf0.*tf1)   ;
        M2=-(tf0.*tf0)  ;
        
       
       tf0(abs(M0) < 1e-10 ) = 0;
       tfc0(cidx, :, tidx) = tf0(1:fLen) ;
     

       omega1 =imag(M1./M0/(2*pi));
       lambda1 =imag(M2./M0/(2*pi));
       lambda2=chirp+lambda1;
       omega2=zeros(size(omega1));
      
      
        for i=1:length(omega1)
           omega2(i)=(i-1)*del_f+omega1(i);
        end
             
       lambda(cidx, :, tidx)=lambda2;
       omega(cidx,:,tidx)=omega2;  
       
    


       lambda(cidx, :, tidx)=lambda2;
       omega(cidx,:,tidx)=omega2; 

    end
  
end
fprintf('\n') ;
end







