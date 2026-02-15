
%%%%%%%%%%%%%%% Example2 %%%%%%%%%%%%%%

% %July, 2025, by Shuixin Li

  %% this code for test the phase  is linear chirp signal in frequency domain
   
   clear;      close all;  clc;      
    Hz = 1024;
    N = 0.5*1024;
    t = (0:N-1)/Hz;
    f = (0:round(N/2))*Hz/N;
    %Mode1
   A1 =exp(-0.00002*(f-256).^2) ;
    %A1 = ones(1,length(f));
    Phi1 = 0.0003*f.*f+0.1*(f);
    GD1 = 0.0006*f+0.1;
    GDD1 = 0.0006*ones(size(GD1));
    X1 = A1.*exp(-1i*2*pi*Phi1);
    X1(end) = -A1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);

    %Mode2
    A2 =exp(-0.00003*(f-256).^2) ;
    %A2 = ones(1,length(f));
    Phi2 = -0.0002*f.*f+0.356*f;
    GD2 = -0.0004*f+0.356;

    GDD2 =  -0.0004*ones(size(GD1)) ;
    X2 = A2.*exp(-1i*2*pi*Phi2);
    X2(end) = -A2(end);
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
    %Test Signal
 
% 
y = (y1+y2); % 保存原始无噪声信号



gs_del=0.1;
gs=10:gs_del:30; arange=0.0016;

SNR_levels = 0:1:50;
num_levels = length(SNR_levels);

% noisy_signals = zeros(length(SNR_levels), length(y));
% for i = 1:length(SNR_levels)
%     noisy_signals(i, :) = awgn(y, SNR_levels(i), 'measured');  % 正确：赋值给第i行
% end
% save('noise_signals.mat');

%num_levels =51; SNR_levels = 0:1:50;
load('noise_chirpsignal.mat')
Sigma_chirp= zeros(1, num_levels);

% Sigma_sinc=sigma;
%save('Sigma_chirp.mat',"Sigma_chirp");
%load('Sigma_chirp.mat')

%%%%%%%%%%%%%%% Example4 %%%%%%%%%%%%%%

l=2.5;


%t^2g,for sigma2 findsigma for FCT
for i = 1:num_levels

 y =noise_chirpsignal (i, :);
[opt_gs1,entro1]=find_FCT_sigma(y,Hz,gs,arange,l);
 Sigma_chirp(i)=  opt_gs1;

end

