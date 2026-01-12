   % %July, 2025, by Shuixin Li
   clear;      close all;  clc;      
    %% Test Signal
    %Parameters
       Hz = 1024;
    N = 0.5*1024;
    t = (0:N-1)/Hz;
    f = (0:round(N/2))*Hz/N;
    %Mode1
    A1 =exp(-0.00025*(f-256).^2) ;
    %A1 = ones(1,length(f));
    Phi1 =-51.2/pi*sin(pi*f./256)+0.25*(f);
    GD1 = -0.2*cos(pi*f./256)+0.25;
    GDD1 = pi/1280*sin(pi*f./256);
    X1 = A1.*exp(-1i*2*pi*Phi1);
    X1(end) = -A1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);

    %Mode2
    A2 =exp(-0.00032*(f-256).^2) ;
    %A2 = ones(1,length(f));
    Phi2 =51.2/pi*sin(pi*f./256)+0.25*f;
    GD2 = 0.2*cos(pi*f./256)+0.25;

    GDD2 = -pi/1280*sin(pi*f./256) ;
    X2 = A2.*exp(-1i*2*pi*Phi2);
    X2(end) = -A2(end);
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
    %Test Signal
    y =(y1+y2)';
 



sigma= 0.007;


chrrange=2500;

tic
[tfc0, lambda1,omega1,tfrtic, tcrtic] =CT_2nd_GPU (y, Hz, sigma,chrrange);                                                        
toc


[tfrsq1] = Msqueeze_CT(y, tfc0, lambda1, omega1, tfrtic,tcrtic, 1);
toc



% tic
% [tfc0_gpu] = CT_gpu(y, Hz, sigma,chrrange);
% toc




%tfc0_chirp = squeeze(tfc0(257, :, :)-tfc0_gpu(257, :, :));  % 得到 L×N 的矩阵
% surf(t, f, abs(tfc0_chirp), 'EdgeColor', 'none');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% zlabel('Amplitude');
% title('Chirplet Transform Amplitude');
% 





% tic
% [tfc0,tfc1, tfc2, lambda1,omega1,tfrtic, tcrtic] = CT_gpu_masked(y, Hz, sigma,chrrange);
% toc

% aa=max(abs(tfc000-tfc0));



del_f=tfrtic(2)-tfrtic(1);
del_c=tcrtic(2)-tcrtic(1);
%[tfrsq1] = Msqueeze_CT(y, tfc0, lambda1, omega1, tfrtic,tcrtic, 1);
tfrsq2 = permute(tfrsq1, [1, 3, 2]);

 [GDs,GDDs]=ridge_3d(tfrsq2,2,30,30,30);

 for i=1:length(GDs)    
       gd01 = t(GDs(:,1));
       gd02 =t(GDs(:,2));
       gdd01=1./tcrtic(GDDs(:,1));
       gdd02=1./tcrtic(GDDs(:,2));         
 end


%  group=zeros(length(f),2);
% for  fidx=1:length(f)
%     group(fidx,1)=R(GDDs(fidx,1),GDs(fidx,1),fidx);
%     group(fidx,2)=R(GDDs(fidx,2),GDs(fidx,2),fidx);
% end
% aaa=real(group);bbb=imag(group);
% 
% plot(f,aaa(:,1))
% plot(f,tcrtic(GDDs(:,1)))

if sum((gd01-GD1).^2)<sum((gd01-GD2).^2)
    gd1=gd01;gd2=gd02;
else 
    gd1=gd02;gd2=gd01;
end

if sum((gdd01-GDD1).^2)<sum((gdd01-GDD2).^2)
    gdd1=gdd01;gdd2=gdd02;
else 
   gdd1=gdd02;gdd2=gdd01;
end

BBB=abs(squeeze(tfc0(:,129,:)));
figure
imageSQ(t, tcrtic, BBB);
shading interp; 
axis xy; 
xlabel('Time (s)', 'FontSize', 20);      
ylabel('Chirprate (Hz/s)', 'FontSize', 20); 


%% Plot of the extracted group decay
figure
e_pdph1 = plot(gd1, f, 'r:', 'linewidth', 2);  % 时间为横轴，频率为纵轴
hold on
pdph1 = plot(GD1, f, 'b-', 'linewidth', 1);
legend([e_pdph1, pdph1], 'Estimated GD1', 'Ground truth', 'Location', 'northwest');
set(legend, 'FontName', 'Times New Roman')
xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
set(gca, 'FontSize', 20)
axis xy;  
%axis([t(1) t(end) f(1) f(end)]);  

figure
e_pdph2 = plot(gd2, f, 'r:', 'linewidth', 2);  % 时间为横轴，频率为纵轴
hold on
pdph2 = plot(GD2, f, 'b-', 'linewidth', 1);
legend([e_pdph2, pdph2], 'Estimated GD2', 'Ground truth', 'Location', 'northeast');
set(legend, 'FontName', 'Times New Roman')
xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
set(gca, 'FontSize', 20)
axis xy;  
%axis([t(1) t(end) f(1) f(end)]); 


figure
e_pdph1 = plot(gdd1, f, 'r:', 'linewidth', 2);  % 时间为横轴，频率为纵轴
hold on
pdph1 = plot(GDD1, f, 'b-', 'linewidth', 1);
legend([e_pdph1, pdph1], 'Estimated GDD1', 'Ground truth', 'Location', 'northwest');
set(legend, 'FontName', 'Times New Roman')
xlabel('GDD(s/Hz)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
set(gca, 'FontSize', 20)
axis xy;  
%axis([0.0001 0.0008 f(1) f(end)]); 


figure
e_pdph2 = plot(gdd2, f, 'r:', 'linewidth', 2);  % 时间为横轴，频率为纵轴
hold on
pdph2 = plot(GDD2, f, 'b-', 'linewidth', 1);
legend([e_pdph2, pdph2], 'Estimated GDD2', 'Ground truth', 'Location', 'northeast');
set(legend, 'FontName', 'Times New Roman')
xlabel('GDD (s/Hz)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
set(gca, 'FontSize', 20)
axis xy; 
%axis([-0.0008 0.0002 f(1) f(end)]);



