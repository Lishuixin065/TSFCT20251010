% %%%%%%%%%%%%%%%% Example1 %%%%%%%%%%%%%%
%    clear;      close all;  clc;      
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
%     y =(y1+y2);
%     gs_del=0.1;
%     gs=10:gs_del:40; chrrange=0.0016; %% 出来的结果挺好
%     %l=2.5 gs=24.6-25;  l=3,=32-35
%%%%%%%%%%%%%%%% Example2 %%%%%%%%%%%%%%
  % clear;close all;  clc;      
  % 
  %  Hz = 512;
  %   N = 1*Hz;
  %   t = (0:N-1)/Hz;
  %   f = (0:round(N/2))*Hz/N;
  %   %Mode1
  %   A1 =exp(-0.000025*(f-128).^2) ;
  %   %A1 = ones(1,length(f));
  % 
  %   Phi1 =0.35*f+(0.2/256^2)*(f.^3);
  %   GD1 =0.6/(256)^2*(f).^2 +0.35;
  %   GDD1 = (1.2/256^2)*f ;
  %   X1 = A1.*exp(-1i*2*pi*Phi1);
  %   X1(end) = -A1(end);
  %   Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
  %   y1 = ifft(Y1);
  % 
  %   %Mode2
  %   A2 =exp(-0.000032*(f-128).^2) ;
  %   %A2 = ones(1,length(f));
  %   Phi2 =-0.2/256^2*f.^3+0.65*f;
  %   GD2 = -0.6/(256)^2*(f).^2 +0.65;
  %   GDD2 =-(1.2/256^2) * f;
  %   X2 = A2.*exp(-1i*2*pi*Phi2);
  %   X2(end) = -A2(end);
  %   Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
  %   y2 = ifft(Y2);
  %   %Test Signal
  %   y =(y1+y2);
  %   chrrange=0.012;
  %    gs_del=0.1;
  %    gs=7:gs_del:9;%7.7-7.8 
%%%%%%%%%%%%%%% Example3 %%%%%%%%%%%%%%
   clear;      close all;  clc;      
    Hz = 1024;
    N = 0.5*1024;
    t = (0:N-1)/Hz;
    f = (0:round(N/2))*Hz/N;
    %Mode1
    A1 =exp(-0.00025*(f-256).^2) ;
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
    y =(y1+y2); arange=0.01;
   gs_del=0.1;
  gs=16:gs_del:20; chrrange=0.01; %% 出来的结果挺好
  %l=2.5 gs=17.1; l=3,25

%%%%%%%%%%%%%%% Example4 %%%%%%%%%%%%%%










l=2.5;

%t^2g,for sigma2 findsigma for FCT
[opt_gs1,entro1]=find_FCT_sigma(y,Hz,gs,chrrange,l);
