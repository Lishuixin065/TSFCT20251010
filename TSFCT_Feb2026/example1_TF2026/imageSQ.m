% imageSQ.m
%
% Display time-frequency result of the Synchrosqueezing transform and others
% function imageSQ(t, ytic, M) ;
%
function imageSQ(t, ytic, M) 

fz = 20;

S = size(M);
Q = M(:);
q = quantile(Q, 0.999);
M(find(M>q)) = q;

imagesc(t, ytic, M)
axis xy ;
set(gca, 'fontsize', fz);
ylabel('Frequency (Hz)', 'FontSize', 20)  
xlabel('Time (s)', 'FontSize', 20)       
%ylim([20,60])