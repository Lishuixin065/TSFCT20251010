
tfproj1 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj1,1)
    for j = 1:size(tfproj1,2)
        tfproj1(i,j) = sum((abs(tfrsq1(:,i,j)).^(2)));
    end
end
figure
imageSQ( t0,tfrtic, tfproj1');
  shading interp; 
  xlabel('Time (s)', 'FontSize', 20);      % 横轴为时间
  ylabel('Frequency (Hz)', 'FontSize', 20); % 纵轴为频率
   rectangle('Position',[0.2 200 0.1 100],'EdgeColor','red','Linewidth',1);
  set(gca, 'FontSize', 20)



