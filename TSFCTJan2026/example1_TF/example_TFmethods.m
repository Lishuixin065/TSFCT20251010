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
    y =(y1+y2);

    %Test Signal
    s=20.5;

   %% SST_2nd
 tic
  [Tx1,tfr,f,IF] = SST_2nd(y, Hz,s);
   fprintf('SST_2nd execution time: ');
toc
 
   figure
  imageSQ( t,f, abs(Tx1));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
  axis xy; 

   figure
  imageSQ( t,f, abs(Tx1));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  xlim([0.2 0.3])
  ylim([200,300])
  axis xy;













 %% TSST_2nd 
 tic
   [Tx2,tfc0,t,f,GroupDelay] = TSST_2nd(y,Hz,s);
 fprintf('TSST_2nd execution time: ');
 toc

   figure
  imageSQ( t,f, abs(Tx2'));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
  axis xy;

  figure
  imageSQ( t,f, abs(Tx2'));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  xlim([0.2 0.3])
  ylim([200,300])
  axis xy;




%% TET_2nd 
tic
   [Tx3,tfc0,t,f,GroupDelay] = TET_2nd(y,Hz,s);
   fprintf('TET_2nd execution time: '); 
toc
   figure
  imageSQ( t,f, abs(Tx3'));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
  axis xy;


 figure
  imageSQ( t,f, abs(Tx3'));
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  xlim([0.2 0.3])
  ylim([200,300])
  axis xy;

%% TSFCT 
arange=0.0016;
[tfrsq1, tfc0, tfrtic, tcrtic] = SFCT2_with_squeeze_gpu(y, Hz, s, arange, 1, 0.01);
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for m = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(m,j) = sum((abs(tfrsq1(:,m,j)).^(2)));
    end
end

 figure
  imageSQ( t,f, abs(tfproj)');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  axis xy;

 figure
  imageSQ( t,f, abs(tfproj)');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  xlim([0.2 0.3])
  ylim([200,300])
  axis xy;


  l=2.5;
  entro_SST_2nd=renyi(Tx1,l);
  entro_TSST_2nd=renyi(Tx2,l);
  entro_TET_2nd=renyi(Tx3,l);
   entro_TSFCT=renyi(tfproj,l);
fprintf('Renyi Entropy Results:\n');
fprintf('  SST-2nd:   %.4f\n', entro_SST_2nd);
fprintf('  TSST-2nd:  %.4f\n', entro_TSST_2nd);
fprintf('  TET-2nd:   %.4f\n', entro_TET_2nd);
fprintf('  TET-2nd:   %.4f\n', entro_TSFCT);

