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

%% plot GDs and GDDs
figure
pddph1=plot(GDD1,f,'r-','linewidth',2);
hold on
pddph2=plot(GDD2,f,'b-','linewidth',2);
legend([pddph1,pddph2],'\rho_1^{\prime\prime}(\eta)','\rho_2^{\prime\prime}(\eta)','Location', 'northwest');
ylabel('Frequency (Hz)', 'FontSize', 20);  
xlabel('','FontSize',20);
set(gca,'FontSize',20)
axis xy; axis([-0.001 0.001 f(1) f(end)]); 

figure
pddph1 = plot(GD1, f, 'r-', 'linewidth', 2);  % 时间是横轴，频率是纵轴
hold on
pddph2 = plot(GD2, f, 'b-', 'linewidth', 2);
legend([pddph1, pddph2], '\rho_1^{\prime}(\eta)', '\rho_2^{\prime}(\eta)','Location', 'northwest');
xlabel('Time (s)', 'FontSize', 20);       
ylabel('Frequency (Hz)', 'FontSize', 20);  
set(gca, 'FontSize', 20)
axis xy;axis([t(1) t(end) f(1) f(end)]); 


 arange=0.0016;
 sigma=25;%entro=25; l=2.5
 

%% squeezeing transform

% tic
% [tfrsq1, tfc0, tfrtic, tcrtic] = SFCT2_with_squeeze_gpu(y, Hz, sigma, arange, 1, 0.01);
% toc

% tic
% [tfc0,  tfrtic, tcrtic, T, R] = SFCT2_GPU(y, Hz, sigma, arange);
% [tfrsq1] = Msqueeze_FCT(y, tfc0, R, T, t, tcrtic, 1, 0.001);
% toc

% ridge detection
  [GDs,GDDs]=ridge_3d(tfrsq1,2,15,15,15);

 for i=1:length(GDs)    
       gd01 = t(GDs(:,1));
       gd02 =t(GDs(:,2));
       gdd01=tcrtic(GDDs(:,1));
       gdd02=tcrtic(GDDs(:,2));         
 end


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
axis([t(1) t(end) f(1) f(end)]);  

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
axis([t(1) t(end) f(1) f(end)]); 


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
axis([0.0001 0.0008 f(1) f(end)]); 


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
axis([-0.0008 0.0002 f(1) f(end)]);



%% Ideal reconstruction

% % real 
% for j=1:length(f)
%     diff1=abs(t-GD1(j)); diff2=abs(t-GD2(j));  
%     [minValue1, minIndex1] = min(diff1); [minValue2, minIndex2] = min(diff2);
%     GDs(j,1)=minIndex1; GDs(j,2)=minIndex2;
%     diff3=abs(tcrtic-GDD1(j)); diff4=abs(tcrtic-GDD2(j));     
%     [minValue3, minIndex3] = min(diff3); [minValue4, minIndex4] = min(diff4);
%     GDDs(j,1)=minIndex3; GDDs(j,2)=minIndex4;
% end
% plot(f,t(GDs(:,1))-GD1)
% plot(f,tcrtic(GDDs(:,1))-GDD1)
% [recov1] = frecov_SSO(GD1, GD2, GDD1, GDD2,GDs,GDDs,tfc0,sigma);
% plot(f,real(recov1(1,:)-conj(X1)))
% [recov2] = ffrecov(GD1, GD2, GDD1, GDD2,GDs,GDDs,tfc2,sigma);
% plot(f,real(recov2(:,1)-conj(X1)'))


%% Frequency domian SSO
[frecov,norm,cond] = frecov_SSO(gd01, gd02, gdd01, gdd02, GDs,GDDs,tfc0, sigma);
%plot the frequency reconstruction errors
if sum(abs(X1-frecov(1,:)).^2) <sum(abs(X2-frecov(1,:)).^2)
  figure, plot(f, real(X1-frecov(1,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20);
   ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{x}_1$', 'Interpreter', 'latex')
  figure, plot(f, real(X2-frecov(2,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20);
  ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
legend('Recovery error of $\hat{x}_2$', 'Interpreter', 'latex')
else
 figure, plot(f, real(X2-frecov(1,:)), 'color', 'k');
   xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20);
 ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{x}_2$', 'Interpreter', 'latex')
  figure, plot(f, real(X1-frecov(2,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency(Hz)','FontSize',20);
  ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{x}_1$', 'Interpreter', 'latex')
end




%% Reconstruction to time domain
trecov=zeros(2,N);
           if mod(N,2)==0        
                for i = 0:N/2
                    trecov(:,i+1) = frecov(:,i+1);
                end
                for i = N/2+1:N-1
                    trecov(:,i+1) = conj(frecov(:,N-i+1));
                end
            else
                for i = 0:(N-1)/2
                    trecov(:,i+1) =frecov(:,i+1);
                end
                for i = (N+1)/2:N-1
                    trecov(:,i+1) = conj(frecov(:,N-i+1));
                end
            end
            %重构计算
            trecov = ifft(trecov, [], 2);     
            trecov = real(trecov);%去除计算误差
          figure
          plot(t,y2)
          hold on
          plot(t,trecov(1,:)-y1)
 

%% recover mean square error  

 [MSEGD1, MSEGD2, MSEGDD1, MSEGDD2] = MSE_ridge_new(gd1, gd2, gdd1, gdd2, GD1, GDD1, GD2, GDD2);
fprintf('%.6f  %.6f  %.6f  %.6f\n', MSEGD1, MSEGD2, MSEGDD1, MSEGDD2)
% 0.001265  0.001164  0.000003  0.000007



[errorx1,errorx2] = MSE_recov_error(y1,y2,trecov); %use part of [1/8, 7/8]
disp([errorx1,errorx2])
% 0.0710    0.0778

%% plot the time domian reconstruction errors
if sum(abs(y1-trecov(1,:)).^2) <sum(abs(y2-trecov(1,:)).^2)
  figure, plot(t, real(y1-trecov(1,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20); 
  ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of ${x}_1$', 'Interpreter', 'latex')
  figure, plot(t, real(y2-trecov(2,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
  ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
 legend('Recovery error of ${x}_2$', 'Interpreter', 'latex')
else
 figure, plot(t, real(y2-trecov(1,:)), 'color', 'k');
 xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
    ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of ${x}_2$', 'Interpreter', 'latex')
  figure, plot(t, real(y1-trecov(2,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
   ylabel('Amplitude (a.u.)','FontSize',20);
  set(gca,'FontSize',20)
 legend('Recovery error of ${x}_1$', 'Interpreter', 'latex')
end


figure 
plot(f, norm, 'LineWidth', 2) 
xlabel('Frequency (Hz)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('The infty norm ', 'FontSize', 20,'FontName', 'Times New Roman');  
set(gca, 'FontSize', 20);
xlim([0,512]);
ylim([0.8,2.2]); 
axis xy;
 


figure 
plot(f, cond, 'LineWidth', 2) 
xlabel('Frequency (Hz)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Condition number ', 'FontSize', 20,'FontName', 'Times New Roman');  
set(gca, 'FontSize', 20);
xlim([0,512]);
ylim([0.8,3]); 
axis xy;
