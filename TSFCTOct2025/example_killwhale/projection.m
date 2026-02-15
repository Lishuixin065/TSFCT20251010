
tfproj1 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj1,1)
    for j = 1:size(tfproj1,2)
        tfproj1(i,j) = sum((abs(tfrsq1(:,i,j)).^(2)));
    end
end
figure
imageSQ( t,tfrtic/1000, tfproj1');
shading interp; 
xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (KHz)', 'FontSize', 20); % 纵轴为频率
rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
set(gca, 'FontSize', 20)
axis square; 


tfproj2 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj2,1)
    for j = 1:size(tfproj2,2)
        tfproj2(i,j) = max(abs((tfrsq1(:,i,j)).^(2)));
    end
end
figure
imageSQ( t,tfrtic/1000, tfproj2');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (KHz)', 'FontSize', 20); % 纵轴为频率
  set(gca, 'FontSize', 20)
axis square; 
  % figure
  % imageSQ( t,tfrtic, abs(tfproj2)');
  % shading interp; 
  % xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  % ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
  % set(gca, 'FontSize', 20)
  % xlim([0.2 0.3])
  % ylim([200,300])
  % axis xy;



  t_idx = 1:length(t);

% 创建一个新的图形窗口
figure;
set(gcf, 'Position', [100, 100, 600, 400]); % 设置图形窗口的大小为 600x400 像素

QN = 5;
D = tfrsq1(:,t_idx,:);
thresh = quantile(abs(D(:)),0.995);
D(abs(D) < thresh * (10-QN+1)/10) = thresh * (10-QN)/10;

hold on; % 保持当前图形，以便在同一图上绘制多个散点图
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10);
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3(t(I2 ), tfrtic(I3)/1000 , tcrtic(I1) , 10,[1 1 1]*(jj-1)/8 , 'filled');   
end
shading interp;
view(0,90); % 设置视角为从上往下看

% 设置坐标轴标签和字体大小
xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
ylabel('Frequency (KHz)', 'FontSize', 20); % 纵轴为频率
set(gca, 'FontSize', 20); % 设置坐标轴字体大小
ylim([0 16])
% 添加颜色条（如果需要）
%colorbar;

% 保存图形（如果需要）
% saveas(gcf, '一致的图形.png');
