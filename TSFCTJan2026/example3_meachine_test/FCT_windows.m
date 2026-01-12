function [gt] =FCT_windows(gs,ch,type)
      L0=1+1i*2*pi*ch*gs^2;
  if type==0
          gt=@(t) 1./sqrt(L0).*exp(-2*pi^2*gs^2*t.^2./L0);
  elseif type==1   
         gt=@(t) -1i*2*pi*gs^2*t.*(L0).^(-3/2).*exp(-2*pi^2*gs^2*t.^2./L0);
  elseif  type==2     
          gt = @(t) gs^2*((L0).^(-3/2)-((2*pi*gs*t).^2).*(L0).^(-5/2)).*exp(-2*pi^2*gs^2*t.^2./L0);   
 elseif  type==3     
          gt = @(t) -1i*2*pi*gs^4*t.*(3*(L0).^(-5/2)-((2*pi*gs*t).^2).*(L0).^(-7/2)).*exp(-2*pi^2*gs^2*t.^2./L0);  
elseif  type==4     
          gt = @(t) gs^4*(3*(L0).^(-5/2)-6*((2*pi*gs*t).^2).*(L0).^(-7/2)+((2*pi*gs*t).^4).*(L0).^(-9/2)).*exp(-2*pi^2*gs^2*t.^2./L0);        
 
 elseif  type==5     
          gt = @(t) gs^6*((-30*pi*1i*t).*(L0).^(-7/2)+80*(1i*pi^3*t.^3*gs^2).*(L0).^(-9/2)-(32*1i*pi^5*t.^5*gs^4).*(L0).^(-11/2)).*exp(-2*pi^2*gs^2*t.^2./L0);

  elseif  type==6     
          gt = @(t) gs^6*(15*(L0).^(-7/2)-180*(pi^2*t.^2*gs^2).*(L0).^(-9/2)+(240*pi^4*t.^4*gs^4).*(L0).^(-11/2)-((2*pi*t.*gs).^6).*(L0).^(-13/2)).*exp(-2*pi^2*gs^2*t.^2./L0); 
  else
            error('Unknown window type: %s', type);
   end 
end



