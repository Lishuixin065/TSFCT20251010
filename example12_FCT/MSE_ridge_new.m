
function [MSEIFx1, MSEIFx2, MSECRx1, MSECRx2] = MSE_ridge_new(if1, if2, chirp1, chirp2, der_x1, dder_x1, der_x2, dder_x2)
% MSE_ridge Calculates the mean square error (MSE) of the instantaneous frequency and chirp rate.
% Input parameters:
% if1, if2 : Estimated instantaneous frequencies for two signal sets
% chirp1, chirp2 : Estimated chirp rates for two signal sets
% der_x1, dder_x1 : True instantaneous frequency and chirp rate of signal x1
% der_x2, dder_x2 : True instantaneous frequency and chirp rate of signal x2
% Output parameters:
% MSEIFx1 : MSE of the instantaneous frequency for signal x1
% MSEIFx2 : MSE of the instantaneous frequency for signal x2
% MSECRx1 : MSE of the chirp rate for signal x1
% MSECRx2 : MSE of the chirp rate for signal x2

len=length(der_x1);
tt=round(len./8):round(7*len./8);ttlen=length(tt);
l=1;
% Determine which estimated values have the smaller error and calculate the MSE accordingly
if  sum((abs(if1(tt) - der_x1(tt))).^l)< sum((abs(if2(tt) - der_x1(tt))).^l)
    MSEIFx1 = 1/(ttlen)*sum((abs(if1(tt) - der_x1(tt))./der_x1(tt)).^l);
    MSEIFx2 =  1/(ttlen)*sum((abs(if2(tt) - der_x2(tt))./der_x2(tt)).^l);
else
    MSEIFx1 =  1/(ttlen)*sum((abs(if2(tt) - der_x1(tt))./der_x1(tt)).^l);
    MSEIFx2 = 1/(ttlen)* sum((abs(if1(tt) - der_x2(tt))./der_x2(tt)).^l); 
end

if  sum((abs(chirp1(tt) - dder_x1(tt))).^l)<sum((abs(chirp2(tt) - dder_x1(tt))).^l) 
    MSECRx1 = 1/(ttlen)* sum((abs(chirp1(tt) - dder_x1(tt))./der_x1(tt)).^l);
    MSECRx2 = 1/(ttlen)*sum((abs(chirp2(tt) - dder_x2(tt))./der_x2(tt)).^l);
else
    MSECRx1 = 1/(ttlen)* sum((abs(chirp2(tt) - dder_x1(tt))./der_x1(tt)).^l);
    MSECRx2 = 1/(ttlen)* sum((abs(chirp1(tt) - dder_x2(tt))./der_x2(tt)).^l); 
end

