function [tfr1, tfrtic, tcrtic] = CT_t2(x, Hz, sigma,chrrange)


[xrow,xcol] = size(x) ;
	%% check input signals
if (xcol~=1)
    error('x must have only one column');

end
Lh=xrow;N=xrow;
crate=linspace(-chrrange,chrrange,2*round(N/2)+1);

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

tcrtic = crate;% discretization of chirp rate

tLen = length(t(1:length(x))) ;% number of bins of time
fLen =length(tfrtic) ; % number of bins of frequency
cLen = length(crate); % number of bins of chirp rate


	%% run STFT and reassignment rule
tfr1=zeros(cLen, fLen, tLen); % synchrosqueezed chirplet transform

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
        tf0 = zeros(N, 1) ;
        tf0(indices) = x(ti+tau).*(time(Lh+1+tau)').^2.*g(Lh+1+tau).*exp(-pi*1i*chirp.*(time(Lh+1+tau)').^2); % for CT with window g          
        tf0 = fft(tf0)/Hz ; tf0 = tf0(1:fLen) ;     
        tfr1(cidx, :, tidx) = tf0(1:fLen) ; 

   
    end
  
end
fprintf('\n') ;
end







