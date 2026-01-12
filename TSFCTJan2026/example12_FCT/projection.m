
tfproj1 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj1,1)
    for j = 1:size(tfproj1,2)
        tfproj1(i,j) = sum((abs(tfrsq1(:,i,j)).^(2)));
    end
end
figure
imageSQ( t,tfrtic, tfproj1');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
   rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
  set(gca, 'FontSize', 20)



% tfproj2 = zeros(size(tfrsq1,2),size(tfrsq1,3));
% t_idx = 1:length(t);
% for i = 1:size(tfproj2,1)
%     for j = 1:size(tfproj2,2)
%         tfproj2(i,j) = max(abs((tfrsq1(:,i,j)).^(2)));
%     end
% end
% figure
% imageSQ( t,tfrtic, tfproj2');
%   shading interp; 
%   xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
%   ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
%   set(gca, 'FontSize', 20)

  figure
  imageSQ( t,tfrtic, abs(tfproj1)');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
  xlim([0.2 0.3])
  ylim([200,300])
  axis xy;