% 主程序保持不变，只修正STFT函数
clear;
% 1. 指定文件名
filename = '1-s2.0-S0888327017300584-mmc1.mp3';

% 2. 读取音频文件
[data, Fs] = audioread(filename);

% 3. 显示音频信息
info = audioinfo(filename);
fprintf('文件信息:\n');
fprintf('文件名: %s\n', info.Filename);
fprintf('采样率: %d Hz\n', info.SampleRate);
fprintf('总时长: %.2f 秒\n', info.Duration);
fprintf('声道数: %d\n', info.NumChannels);
fprintf('总采样点数: %d\n', info.TotalSamples);

% 4. 处理多声道音频（转换为单声道）
if size(data, 2) == 2
    fprintf('这是立体声音频，正在转换为单声道...\n');
    data = mean(data, 2);
else
    fprintf('这是单声道音频。\n');
end

% 6. 保存数据（可选）
save('audio_data.mat', 'data', 'Fs');
fprintf('音频数据已保存为 audio_data.mat\n');

fprintf('\n转换完成！音频数据已存储在变量 data 中，采样率存储在变量 Fs 中。\n');

% Parameter setting
fprintf('Initialization....\n')

f_start = 1250;
f_end = 2000;
data_band_original = bandpass(data, [f_start, f_end], Fs);

% 修正：对实信号做希尔伯特变换
data_analytic = hilbert(data_band_original);
% 修正：正确进行频移操作
t_original = (0:length(data_band_original)-1)'/Fs;
data = data_analytic .* exp(-1i*2*pi*1250*t_original);

m = 32;
data1 = downsample(data, m); % downsample the signal by 32
fs = Fs/m;

L = 110;
data2 = data1(round(L*fs)+1: round((L+40)*fs));
tt = L + (0:length(data2)-1)/fs;
N = length(data2);

figure;
plot(tt, real(data2)); % 显示实部
xlabel('Time (s)', 'FontSize', 20);
ylabel('Amplitude (a.u.)', 'FontSize', 20);
%title('Signal Waveform ', 'FontSize', 20);


hop=100;


s=1;  %sigma for Frequency-domain Gaussian window




%% STFT 
tic
[SST2,STFT,f] = SST_2nd_hop_gpu(data2,s,hop,fs);
toc
fprintf('STFT... \n ');
ff=f+1300;
figure
imageSQ(tt,ff/1000,abs(STFT));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)



%% SST_2nd


figure
fprintf('SST_2nd... \n ');
imageSQ(tt,ff/1000,abs(SST2));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)



%% TSST
tic
[TSST2,tfc0,t,f,T] = TSST_2nd_hop_gpu(data2,fs,s,hop);
toc
figure
imageSQ(tt,ff/1000,abs(TSST2));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)



%% TET_2nd
tic
[TET2,~,~,~,~] = TET_2nd_hop_gpu(data2,fs,s,hop);
toc
fprintf('TET_2nd... \n ');
figure
imageSQ(tt,ff/1000,abs(TET2));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)




%% TSFCT

arange=0.15;

% tic
% [tfc0,tfc1,tfc2,tfrtic,tcrtic,R,T] = FCT2nd_hop_gpu(data2,fs,s,arange,hop);
% t0=(1:hop:N)/fs;
% [tfrsq1] = Msqueeze_FCT(data2, tfc0, R, T,t0,tcrtic, 1,10^(-8));
% toc

tic
[tfc0,tfrsq1,tfrtic,tcrtic,R,T] = Msqueeze_FCT_gpu(data2,fs,s,arange,hop,1,10^(-8));
toc

tfrsq1(1,:,:)=0;  tfrsq1(end,:,:)=0;                 
tfproj1 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj1,1)
    for j = 1:size(tfproj1,2)
        tfproj1(i,j) = sum((abs(tfrsq1(:,i,j)).^(2)));
    end
end
figure
imageSQ( tt,ff/1000, tfproj1');
shading interp; 
xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (KHz)', 'FontSize', 20); % 纵轴为频率
rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
set(gca, 'FontSize', 20)




 l=2.5;
  entro_SST_2nd=renyi(SST2,l);
  entro_TSST_2nd=renyi(TSST2,l);
  entro_TET_2nd=renyi(TET2,l);
  entro_TSFCT=renyi(tfproj1,l);
fprintf('Renyi Entropy Results:\n');
fprintf('  SST-2nd:   %.4f\n', entro_SST_2nd);
fprintf('  TSST-2nd:  %.4f\n', entro_TSST_2nd);
fprintf('  TET-2nd:   %.4f\n', entro_TET_2nd);
fprintf('  FCT-projection:   %.4f\n', entro_TSFCT);
