function[entro] = renyi(tfr,l)


upp1=sum(sum(abs(tfr).^(2*l)));
dow1=sum(sum(abs(tfr).^2))^(l);


entro=1/(1-l)*log(upp1/dow1);