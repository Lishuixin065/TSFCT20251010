% 将MP3文件转换为MATLAB矩阵
clear all;
% 1. 指定文件名
filename = 'CMST_Wellard_BremerCanyon_2014-2015_orca115.wav';

[data, Fs] = audioread(filename);
figure
tt=(0:length(data)-1)/Fs;
plot(tt, data)
xlim([0,tt(end)])
xlabel('Time (s)', 'FontSize', 20);
ylabel('Amplitude', 'FontSize', 20);
title('Audio Waveform', 'FontSize', 20);      
set(gca, 'FontSize', 20)


m=3;% the m-downsampling
data1 = downsample(data,m); 
fs = Fs/m;% the downsampling  sample frequency
L=0.6;L0=0.5;
hop=20;
data2=data1(L*fs+1:(L+L0)*fs);  %the extracted time segment
t=L+(0:hop:length(data2)-1) / fs;

N=length(data2);

  if mod(N,2)==0
        freqs = (0:hop:round(N/2))*fs/N/1000;
  else
        freqs = (0:hop:round((N-1)/2))*fs/N/1000;
  end



%% Calculating TFRs
fprintf('Calculating TFRs... \n ');



chrrange=3*10^(-4);
s=6.0*10^(1); %sigma for frequency-domain Gaussian window

[tfc0,tfrtic,tcrtic,T,R] = SFCT2_hop(data2,fs,s,chrrange,hop);
[tfrsq1] = Msqueeze0(data2, tfc0, R, T,t,tcrtic, 1,10^(-8));
 

tfrsq2=zeros(size(tfrsq1));
tfrsq2(:,:,100:end)=tfrsq1(:,:,100:end);
