
tfproj1 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj1,1)
    for j = 1:size(tfproj1,2)
        tfproj1(i,j) = sum(abs((tfrsq1(:,i,j)).^(1)));
    end
end
figure
imageSQ(t, tfrtic, tfproj1);
shading interp; 
axis xy; xlabel('time (s)'); ylabel('frequency (Hz)');


tfproj2 = zeros(size(tfrsq1,2),size(tfrsq1,3));
t_idx = 1:length(t);
for i = 1:size(tfproj2,1)
    for j = 1:size(tfproj2,2)
        tfproj2(i,j) = max(abs((tfrsq1(:,i,j)).^(1)));
    end
end
figure
imageSQ(t, tfrtic, tfproj2);
shading interp; 
axis xy; xlabel('time (s)'); ylabel('frequency (Hz)');