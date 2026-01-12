
%  % %July, 2025, by Shuixin Li
   clear;      close all;  clc;  


% %%%%%%%%%%%%%%%% Example1 %%%%%%%%%%%%%%
%     Hz = 1024;
%     N = 0.5*1024;
%     t = (0:N-1)/Hz;
%     f = (0:round(N/2))*Hz/N;
%     %Mode1
%    A1 =exp(-0.00002*(f-256).^2) ;
%     %A1 = ones(1,length(f));
%     Phi1 = 0.0003*f.*f+0.1*(f);
%     GD1 = 0.0006*f+0.1;
%     GDD1 = 0.0006*ones(size(GD1));
%     X1 = A1.*exp(-1i*2*pi*Phi1);
%     X1(end) = -A1(end);
%     Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
%     y1 = ifft(Y1);
% 
%     %Mode2
%     A2 =exp(-0.00003*(f-256).^2) ;
%     %A2 = ones(1,length(f));
%     Phi2 = -0.0002*f.*f+0.356*f;
%     GD2 = -0.0004*f+0.356;
% 
%     GDD2 =  -0.0004*ones(size(GD1)) ;
%     X2 = A2.*exp(-1i*2*pi*Phi2);
%     X2(end) = -A2(end);
%     Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
%     y2 = ifft(Y2);
%     %Test Signal
%     y =(y1+y2)';
% 
% chrrange=8000;
% gs_del=0.001;gs=0.001:gs_del:0.01;

% l=2.5; gs=0.008


%%%%%%%%%%%%%%%% Example2 %%%%%%%%%%%%%%  
    %% Test Signal
    %Parameters
   Hz = 1024;
    N = 0.5*1024;
    t = (0:N-1)/Hz;
    f = (0:round(N/2))*Hz/N;
    %Mode1
    A1 =exp(-0.00025*(f-256).^2);
    %A1 = ones(1,length(f));
    Phi1 =-51.2/pi*sin(pi*f./256)+0.25*(f);
    GD1 = -0.2*cos(pi*f./256)+0.25;
    GDD1 = pi/1280*sin(pi*f./256);
    X1 = A1.*exp(-1i*2*pi*Phi1);
    X1(end) = -A1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);

    %Mode2
    A2 =exp(-0.00032*(f-256).^2) ;
    %A2 = ones(1,length(f));
    Phi2 =51.2/pi*sin(pi*f./256)+0.25*f;
    GD2 = 0.2*cos(pi*f./256)+0.25;

    GDD2 = -pi/1280*sin(pi*f./256) ;
    X2 = A2.*exp(-1i*2*pi*Phi2);
    X2(end) = -A2(end);
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
    %Test Signal
    y =(y1+y2)';
 chrrange=2500;
 gs_del=0.001;gs=0.001:gs_del:0.02;
 l=2.5; %gs=0.007



leng=length(gs); entro=zeros(1,leng); 
for j=1:leng

[tfc, tfrtic, tcrtic] = CT(y, Hz, gs(j),chrrange);
%[tfc, tfrtic, tcrtic] = CT_t2(y, Hz, gs(j),chrrange);
upp1=sum(sum(sum(abs(tfc).^(2*l))));
dow1=sum(sum(sum(abs(tfc).^2)))^(l);
entro(j)=1/(1-l)*log(upp1/dow1);
%disp("j"), 
results = [j, gs(j), entro(j)];
disp(results);
end 

[gs_min1, min_ind1]=min(entro);

opt_gs1=gs(min_ind1);
disp(opt_gs1);


