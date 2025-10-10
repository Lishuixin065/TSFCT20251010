function [recov,A0,B0] = frecov_SSO(gd1, gd2, ahirp1, ahirp2, GDs,GDDs,tfc0, sigma)

% Reconstructs signal components using FGSSO
% in the time-frequency-GDD domain.
%
% Inputs:
%   gd1, gd2    : Group delay estimates for components 
%   ahirp1, ahirp2 : Instantaneous GDDs estimates for components 
%   GDs         : Real Group delays  for components
%   GDDs        : Real Group delay dispersion for components
%   tfc0        : 3D time-frequency-group delay dispersion represention  
%   sigma       : Gaussian window width parameter (scalar > 0)
%
% Output:
%   recov       : Recovered frequency-domain signal components 
%    A0         : The infty norm of inver coefficent matrix 




  A0=zeros(1,length(gd1));
  B0=zeros(1,length(gd1));
  recov = zeros(2,length(gd1)); 
  scale = sigma;
  
   for char = 1:length(gd1)
    tmp = zeros(2,2);
    tmp(1,1) = Gaussian_FAW(0,0,scale);
    tmp(1,2) = Gaussian_FAW(ahirp2(char)-ahirp1(char),gd2(char)-gd1(char),scale);   
    tmp(2,1) =Gaussian_FAW(ahirp1(char)-ahirp2(char),gd1(char)-gd2(char),scale);
    tmp(2,2) =Gaussian_FAW(0,0,scale) ;  
    xtmp = [tfc0(GDDs(char,1),GDs(char,1),char); tfc0(GDDs(char,2),GDs(char,2),char)];   
    recov(:,char) = pinv(tmp + 1*eps) * xtmp;
    A0(1,char)=norm(pinv(tmp + 1*eps), inf);
    B0(1,char)=cond(tmp );
   end



% 
% if nargout == 1
%     varargout{1} = recov;
% elseif nargout == 2
%     varargout{1} = recov;
%     varargout{2} = A0;
% else
%     error('不支持的输出参数个数: %d', nargout);
% end




%% AW
function [FAg] = Gaussian_FAW(GDD,GD, sigma) 
   L=1+1i*2*pi*sigma^2*GDD;
   FAg = L^(-0.5)*exp(-2*pi^2*sigma^2*GD^2./L);   
end

end