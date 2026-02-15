 %% 3d plot of chirplet transform with g_0

t_idx = 1:length(t);
t_show = t(t_idx);
figure
QN = 5 ;
D = tfrsq1(:,:,t_idx);


thresh = quantile(abs(D(:)),0.99995);
 D(abs(D) < thresh * (10-QN+1)/10) = thresh * (10-QN)/10; 

for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);

    scatter3((I3 )/Hz, tfrtic(I2) , tcrtic(I1) , 10,[1 1 1]*(jj-1)/8 , 'filled'); 

    hold on
end
shading interp;

view(45,30)

xlabel('Time (s)' ,'FontSize', 20);
ylabel('Frequency (Hz)', 'FontSize', 20);
zlabel( 'Chirprate(Hz/s)','FontSize', 20);
set(gca, 'Box', 'on'); 
title('3D Plot of SCT ','FontSize', 20); 







