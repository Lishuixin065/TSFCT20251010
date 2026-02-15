%This script visualizes the FCT/HTSFCT transforms in 3D space

t_idx = 1:length(t);


figure

QN = 5;
D = tfc2(:,t_idx,:);
thresh = quantile(abs(D(:)),0.99995);
D(abs(D) < thresh * (10-QN+1)/10) = thresh * (10-QN)/10;

for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10);
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3(t(I2), tfrtic(I3)/1000, tcrtic(I1), 10, [1 1 1]*(jj-1)/8, 'filled');   
    hold on
end

shading interp;
view(45,30)
xlabel('Time (s)' ,'FontSize', 20);
ylabel('Frequency (Hz)', 'FontSize', 20);
zlabel('GDD (s/Hz)', 'FontSize', 20);
set(gca, 'Box', 'on'); 
grid on;    
title('3D Plot of TSFCT ', 'FontSize', 20); 
