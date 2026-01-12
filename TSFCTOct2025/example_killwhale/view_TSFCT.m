tfproj2 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj2,1)
    for j = 1:size(tfproj2,2)
        tfproj2(i,j) = max(abs((tfrsq1(:,i,j)).^(2)));
    end
end

% 首先生成参考图并获取其属性
figure
imageSQ(t, tfrtic/1000, tfproj2');
shading interp; 
xlabel('Time (s)', 'FontSize', 20);
ylabel('Frequency (KHz)', 'FontSize', 20);
set(gca, 'FontSize', 20)

% 获取参考图的位置和坐标轴范围
ref_pos = get(gcf, 'Position');
ref_xlim = xlim;
ref_ylim = ylim;
close(gcf); % 关闭参考图

% 现在生成3D散点图，并设置为与参考图一致
figure('Position', ref_pos) % 使用相同的位置和大小
QN = 5;
D = tfrsq1(:,t_idx,:);
thresh = quantile(abs(D(:)),0.995);
D(abs(D) < thresh * (10-QN+1)/10) = thresh * (10-QN)/10;

for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10);
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3(t(I2), tfrtic(I3)/1000, tcrtic(I1), 10, [1 1 1]*(jj-1)/8, 'filled');   
    hold on
end
shading interp;
view(0,90)

% 设置与参考图相同的坐标轴范围
xlim(ref_xlim); % 时间范围一致
ylim(ref_ylim); % 频率范围一致

% 设置相同的标签和字体
xlabel('Time (s)', 'FontSize', 20);
ylabel('Frequency (KHz)', 'FontSize', 20);
zlabel('GDD', 'FontSize', 20);
set(gca, 'FontSize', 20)

% 移除边框
set(gca, 'Box', 'off', 'TickDir', 'out');
set(gcf, 'Color', 'w');
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);
set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off');
