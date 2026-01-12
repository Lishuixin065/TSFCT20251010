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
SNR_levels = 0:1:50;
num_levels = length(SNR_levels);

% noise_chirpsignal = zeros(length(SNR_levels), length(y));
% for i = 1:length(SNR_levels)
%    noise_chirpsignal(i, :) = awgn(y, SNR_levels(i), 'measured');  % 正确：赋值给第i行
% end
% 
%  save('noise_chirpsignal.mat', 'noise_chirpsignal');




%num_levels =51; SNR_levels = 0:1:50;
load('noise_chirpsignal.mat')
load('STFT_sigma.mat')
load('Sigma_chirp.mat')
 Renyi_SST2  = zeros(1, num_levels);
 Renyi_TSST2= zeros(1, num_levels);
 Renyi_TET2= zeros(1, num_levels);
 Renyi_TSFCT = zeros(1, num_levels);

    %Test Signal


   %% SST_2nd
for i = 1:num_levels
    fprintf('Processing SNR = %d dB...\n', SNR_levels(i));

    y = noise_chirpsignal(i, :);
       s=STFT_sigma(i);
  [Tx1,tfr,~,IF] = SST_2nd(y, Hz,s);
  [Tx2,~,~,~,~] = TSST_2nd(y,Hz,s);
  [Tx3,~,t,f,GroupDelay] = TET_2nd(y,Hz,s);

arange=0.0016;
sigma=Sigma_chirp(i);
[tfrsq1, tfc0, tfrtic, tcrtic] = SFCT2_with_squeeze_gpu(y, Hz, sigma, arange, 1, 0.01);
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for m = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(m,j) = sum((abs(tfrsq1(:,m,j)).^(2)));
    end
end

  l=2.5;
  entro_SST_2nd=renyi(Tx1,l);
  entro_TSST_2nd=renyi(Tx2,l);
  entro_TET_2nd=renyi(Tx3,l);
  entro_TSFCT_2nd=renyi(tfproj,l);

  Renyi_SST2(i)  =entro_SST_2nd;
 Renyi_TSST2(i)= entro_TSST_2nd;
 Renyi_TET2(i)= entro_TET_2nd;
 Renyi_TSFCT(i) = entro_TSFCT_2nd;



end

figure;
plot(SNR_levels,Renyi_TSFCT, 'm-d', 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'TSFCT projection');
hold on;
plot(SNR_levels, Renyi_TSST2, 'b-s', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'TSST-2');
plot(SNR_levels, Renyi_TET2, 'g-^', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'TET-2');
plot(SNR_levels, Renyi_SST2, 'r-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'SST-2');
xlim([0, 50]);
xlabel('Input-SNR (dB)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Renyi entroy (bits)', 'FontSize', 20,'FontName', 'Times New Roman');  
legend('show', 'Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 20);