
function [errorx1,errorx2] = MSE_recov_error(x1,x2,recov)
% MSE_recov Calculates the mean square error of signal recovery for two signals.
% Input parameters:
% x1 : Original signal 1
% x2 : Original signal 2
% reco : Recovered signal (matrix with two columns, one for each original signal)
% 
% Output parameters:
% errorx1 : Mean square error of the recovery for signal 1
% errorx2 : Mean square error of the recovery for signal 2


tn=length(x1);
time=ceil(tn/8):floor(tn*7/8);
tlen=length(time);
x1=x1(1,time);
x2=x2(1,time);
recov=recov(:,time);
plot(time,x2-recov(1,:))

l=2;
if sum(abs(x1(:)-recov(1,:)),"all")<sum(abs(x2(:)-recov(1,:)),"all")
  errorx1=sqrt(1/(tlen)*sum((abs(x1-recov(1,:))).^l)/sum(abs(x1.^l)));
  errorx2=sqrt(1/(tlen)*sum((abs(x2-recov(2,:))).^l)/sum(abs(x2.^l)));
else
  errorx1=sqrt(1/(tlen)*sum((abs(x1-recov(2,:))).^l)/sum(abs(x1.^l)));
  errorx2=sqrt(1/(tlen)*sum((abs(x2-recov(1,:))).^l)/sum(abs(x2.^l)));
end






