% 将MP3文件转换为MATLAB矩阵
clear;
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

plot(fft(data2))

N=length(data2);

  if mod(N,2)==0
        freqs = (0:hop:round(N/2))*fs/N/1000;
  else
        freqs = (0:hop:round((N-1)/2))*fs/N/1000;
  end



%% Calculating TFRs
fprintf('Calculating TFRs... \n ');

sigma=2*10^(-3); %sigma for Gaussian window
% STFT 
[STFT,SST2] = SST_2nd_hop(data2,fs,sigma,hop);
fprintf('STFT... \n ');
figure
imageSQ(t,freqs,abs(STFT));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)

% SST_2nd
figure
fprintf('SST_2nd... \n ');
imageSQ(t,freqs,abs(SST2));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)

% TSST_2nd
s=6*10^1;  %sigma for Frequency-domain Gaussian window
[TSST2,~,~,f,T] = TSST_2nd_hop(data2,fs,s,hop);
fprintf('TSST_2nd... \n ');
figure
imageSQ(t,freqs,abs(TSST2').^2);
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)





% TET_2nd
[TET2,~,~,~,~] = TET_2nd_hop(data2,fs,s,hop);
fprintf('TET_2nd... \n ');
figure
imageSQ(t,freqs,abs(TET2'));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)


