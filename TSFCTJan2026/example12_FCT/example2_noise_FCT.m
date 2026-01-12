
% Generate noisy signals with varying SNR levels
clear; close all; clc;
 %Parameters
    Hz = 1024;
    N = 0.5*1024;
    t = (0:N-1)/Hz;
    f = (0:round(N/2))*Hz/N;


    %Mode1
    A1 =exp(-0.00032*(f-256).^2) ;
    %A1 = ones(1,length(f));
    Phi1 =-51.2/pi*sin(pi*f./256)+0.25*(f);
    GD1 = -0.2*cos(pi*f./256)+0.25;
    GDD1 = pi/1280*sin(pi*f./256);
    X1 = A1.*exp(-1i*2*pi*Phi1);
    X1(end) = -A1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);

    %Mode2
    A2 =exp(-0.00025*(f-256).^2) ;
    %A2 = ones(1,length(f));
    Phi2 =51.2/pi*sin(pi*f./256)+0.25*f;
    GD2 = 0.2*cos(pi*f./256)+0.25;

    GDD2 = -pi/1280*sin(pi*f./256) ;
    X2 = A2.*exp(-1i*2*pi*Phi2);
    X2(end) = -A2(end);
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
% 
y = (y1+y2); % 保存原始无噪声信号
SNR_levels = 0:1:50;
num_levels = length(SNR_levels);

% noise_sincsignal = zeros(length(SNR_levels), length(y));
% for i = 1:length(SNR_levels)
%    noise_sincsignal (i, :) = awgn(y, SNR_levels(i), 'measured');  % 正确：赋值给第i行
% end
% save('noise_sincsignal.mat','noise_sincsignal.mat');
% 

%num_levels =51; SNR_levels = 0:1:50;
load('noise_sincsignal.mat')


load('Sigma_sinc.mat')
 RQFnoise_IF  = zeros(1, num_levels);
 RQFnoise_CR= zeros(1, num_levels);
 RQFnoise_sigma= zeros(1, num_levels);
 RQFnoise_sigma0 = zeros(1, num_levels);


 arange=0.01;

for i = 1:num_levels
    fprintf('Processing SNR = %d dB...\n', SNR_levels(i));

    y =noise_sincsignal(i, :);
    sigma=Sigma_sinc(i);
[tfrsq1, tfc0, tfrtic, tcrtic] = SFCT2_with_squeeze_gpu(y, Hz, sigma, arange, 1, 0.01);

[GDs,GDDs]=ridge_3d(tfrsq1,2,20,20,20);

gd01_all = t(GDs(:, 1));    % 第一列所有数据
gd02_all = t(GDs(:, 2));    % 第二列所有数据
gdd01_all = tcrtic(GDDs(:, 1));  % GDDs第一列所有数据
gdd02_all = tcrtic(GDDs(:, 2));  % GDDs第二列所有数据
window_size = 15;  % 
gd01 = smoothdata(gd01_all, 'gaussian', window_size);
gd02 = smoothdata(gd02_all, 'gaussian', window_size);
gdd01 = smoothdata(gdd01_all, 'gaussian', window_size);
gdd02 = smoothdata(gdd02_all, 'gaussian', window_size);       


%% recov_with_sigma
[frecov1] = frecov_SSO(gd01, gd02, gdd01, gdd02, GDs,GDDs,tfc0, sigma);

 trecov1=zeros(2,N);
           if mod(N,2)==0        
                for ii = 0:N/2
                    trecov1(:,ii+1) = frecov1(:,ii+1);
                end
                for ii = N/2+1:N-1
                    trecov1(:,ii+1) = conj(frecov1(:,N-ii+1));
                end
            else
                for ii = 0:(N-1)/2
                    trecov1(:,ii+1) =frecov1(:,ii+1);
                end
                for ii = (N+1)/2:N-1
                    trecov1(:,ii+1) = conj(frecov1(:,N-ii+1));
                end
            end
            %重构计算
            trecov1 = ifft(trecov1, [], 2);     
            trecov1 = real(trecov1);%去除计算误差
       


%% recov_with_sigma0
sigma0=sigma/3;
[tfc] = FCT2_gpu(y,Hz,sigma0,arange);
[frecov2] = frecov_SSO(gd01, gd02, gdd01, gdd02, GDs,GDDs,tfc, sigma0);
 trecov2=zeros(2,N);
           if mod(N,2)==0        
                for ii = 0:N/2
                    trecov2(:,ii+1) = frecov2(:,ii+1);
                end
                for ii = N/2+1:N-1
                    trecov2(:,ii+1) = conj(frecov2(:,N-ii+1));
                end
            else
                for ii = 0:(N-1)/2
                    trecov2(:,ii+1) =frecov2(:,ii+1);
                end
                for ii = (N+1)/2:N-1
                    trecov2(:,ii+1) = conj(frecov2(:,N-ii+1));
                end
            end
            %重构计算
            trecov2 = ifft(trecov2, [], 2);     
            trecov2 = real(trecov2);%去除计算误差
       

 [RQF_IF, RQF_CR] = RQF_ridge(gd01, gd02, gdd01, gdd02, GD1, GDD1, GD2, GDD2);    

 [RQFr1] = RQF_recov(y1, y2, trecov1);
 [RQFr2] = RQF_recov(y1, y2, trecov2);

 RQFnoise_IF(i)  = RQF_IF;
 RQFnoise_CR(i)= RQF_CR;
 RQFnoise_sigma(i)=RQFr1;
 RQFnoise_sigma0(i) =RQFr2;
end



figure;
plot(SNR_levels, RQFnoise_IF, 'r-o', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'GD');
hold on;
plot(SNR_levels, RQFnoise_CR, 'b-s', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'GDD');
xlim([0, 50]);
xlabel('Input-SNR (dB)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('RQF (dB)', 'FontSize', 20,'FontName', 'Times New Roman'); 
legend('show', 'Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 20);



figure;
plot(SNR_levels, RQFnoise_sigma, 'g-^', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'FGSSO with \sigma');
hold on
plot(SNR_levels,RQFnoise_sigma0, 'm-d', 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'FGSSO with \sigma_0');
xlim([0, 50]);
xlabel('Input-SNR (dB)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('RQF (dB)', 'FontSize', 20,'FontName', 'Times New Roman'); 
legend('show', 'Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 20);
