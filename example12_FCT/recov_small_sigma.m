%% Frequency recov
sigma0=sigma;
[tfc00] = FCT(y,Hz,sigma0,arange);
[frecov] = frecov_SSO(gd01, gd02, gdd01, gdd02, GDs,GDDs,tfc00, sigma0);
%plot the frequency reconstruction errors
if sum(abs(X1-frecov(1,:)).^2) <sum(abs(X2-frecov(1,:)).^2)
  figure, plot(f, real(X1-frecov(1,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20); 
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{y}_1$', 'Interpreter', 'latex')
  figure, plot(f, real(X2-frecov(2,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20);
  set(gca,'FontSize',20)
legend('Recovery error of $\hat{y}_2$', 'Interpreter', 'latex')
else
 figure, plot(f, real(X2-frecov(1,:)), 'color', 'k');
   xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency (Hz)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{y}_2$', 'Interpreter', 'latex')
  figure, plot(f, real(X1-frecov(2,:)), 'color', 'k');
  xlim([f(1) f(end)]);
  ylim([-1 1])
  xlabel('Frequency(Hz)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of $\hat{y}_1$', 'Interpreter', 'latex')
end





%% time-recov
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
         
 

%% recover mean square error
 [MSEGD1, MSEGD2, MSEGDD1, MSEGDD2] = MSE_ridge_new(gd1, gd2, gdd1, gdd2, GD1, GDD1, GD2, GDD2);
disp([MSEGD1, MSEGD2, MSEGDD1, MSEGDD2])
% 2nd 0.0023    0.0012    0.0003    0.0002
% 3rd 0.0021    0.0004    0.0001    0.0001   0.0143    0.0147
% 4th 0.0027    0.0003    0.0001    0.0000   0.0143    0.0147

[errorx1,errorx2] = MSE_recov_error(y1,y2,trecov); %use part of [1/8, 7/8]
disp([errorx1,errorx2])


%% plot error
if sum(abs(y1-trecov(1,:)).^2) <sum(abs(y2-trecov(1,:)).^2)
  figure, plot(t, real(y1-trecov(1,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20); 
  set(gca,'FontSize',20)
  legend('Recovery error of ${y}_1$', 'Interpreter', 'latex')
  figure, plot(t, real(y2-trecov(2,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
  set(gca,'FontSize',20)
 legend('Recovery error of ${y}_2$', 'Interpreter', 'latex')
else
 figure, plot(t, real(y2-trecov(1,:)), 'color', 'k');
 xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
  set(gca,'FontSize',20)
  legend('Recovery error of ${y}_2$', 'Interpreter', 'latex')
  figure, plot(t, real(y1-trecov(2,:)), 'color', 'k');
  xlim([t(1) t(end)]);
  ylim([-0.1 0.1])
  xlabel('Time (s)','FontSize',20);
  set(gca,'FontSize',20)
 legend('Recovery error of ${y}_1$', 'Interpreter', 'latex')
end


