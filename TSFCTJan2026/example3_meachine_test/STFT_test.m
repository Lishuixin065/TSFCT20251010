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

%[tfc0,tfc1,tfc2, t, f]= STFT_hop_with(data2,fs,s,hop);
[tfc0,tfc1,tfc2,t,f] = STFT_modify(data2,fs,s,hop);
ff=f+1300;
figure
imageSQ(t,ff,abs(tfc0));
shading interp;
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (KHz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)