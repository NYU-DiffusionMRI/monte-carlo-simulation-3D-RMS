function X = vgfit(D,t,n)

D = D(:); t = t(:); n = n(:);

Nt = 1e1;
Xt = zeros(Nt,1);
err = zeros(Nt,1);
rng(0);
for i = 1:Nt
    x0 = 10/Nt*i;
    options = optimoptions(@lsqnonlin,'Display','off','Algorithm','Levenberg-Marquardt');
    Xi = lsqnonlin(@(x)costfunction(D,t,n,x),x0,[],[],options);
    Xt(i) = Xi;
    err(i) = sum((D-vgmodel(t,n,Xi)).^2);
end

[~,I] = min(err);
X = Xt(I);

end

function [J, dJ] = costfunction(D,t,n,X)

nt = numel(t);
[DX, dD] = vgmodel(t,n,X);
J = 1/2/nt*sum((DX-D).^2);
dJ = 1/nt*sum((DX-D).*dD);

end

function [D, dD] = vgmodel(t,n,r)

D0 = 2;
td = r^2/D0;
bardelta = n/td; dbardelta = -2*n*D0 / r^3;
barDelta = t/td; dbarDelta = -2*t*D0 / r^3;

N=15; 
b = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644 24.3113   27.4571   30.6019 ...
     33.7462   36.8900   40.0334   43.1766   46.3196 49.4624   52.6050   55.7476   58.8900   62.0323];

s = 0; ds = 0;
for k=1:N
   s = s + (2/(b(k)^6*(b(k)^2-1)))*(-2 + 2*b(k)^2*bardelta + ...
        2*(exp(-b(k)^2*bardelta)+exp(-b(k)^2*barDelta)) - ...
        exp(-b(k)^2*(bardelta+barDelta)) - exp(-b(k)^2*(barDelta-bardelta))); 
   ds = ds +  (2/(b(k)^6*(b(k)^2-1)))*( ...
        2*b(k)^2*dbardelta + ...
        2*exp(-b(k)^2*bardelta) .* (- b(k)^2*dbardelta) + ...
        2*exp(-b(k)^2*barDelta) .* (- b(k)^2*dbarDelta) - ...
        exp(-b(k)^2*(bardelta+barDelta)) .* (- b(k)^2*(dbardelta + dbarDelta)) - ... 
        exp(-b(k)^2*(barDelta-bardelta)) .* (- b(k)^2*(dbarDelta - dbardelta)));           
end
s = s.*D0.*td^3;
D = s./(t-n/3)./n.^2;

ds = ds.*D0.*td^3 + 6*s.*r^5 / D0^2;
dD = ds./(t-n/3)./n.^2;

end
