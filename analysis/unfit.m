function X = unfit(D,t,n,Xrange)

D = D(:); t = t(:); n = n(:);

Nt = 1e3;
Xt = zeros(2,Nt);
err = zeros(Nt,1);
rng(0);
x0t = [Xrange(1)*rand(Nt,1) Xrange(2)*rand(Nt,1)];
parfor i = 1:Nt
%     x0 = [Xrange(1)*rand; Xrange(2)*rand];
    x0 = x0t(i,:);
    options = optimoptions(@lsqnonlin,'Display','off','Algorithm','Levenberg-Marquardt');
    Xi = lsqnonlin(@(x)costfunction(D,t,n,x),x0,[],[],options);
    Xt(:,i) = Xi;
    err(i) = sum((D-unmodel(t,n,Xi)).^2);
end

[~,I] = min(err);
X = Xt(:,I);

end

function [J, dJ] = costfunction(D,t,n,X)

nt = numel(t);
[DX, dD] = unmodel(t,n,X);
J = 1/2/nt*sum((DX-D).^2);

dJ = zeros(2,1);
dJ(1) = 1/nt*sum((DX-D).*dD(:,1));
dJ(2) = 1/nt*sum((DX-D).*dD(:,2));

end

function [D, dD] = unmodel(t,n,X)

w = X(1); td = X(2);
bardelta = n/td; dbardelta = -n / td^3;
barDelta = t/td; dbarDelta = -t / td^3;

s = 1/4 * ( 2*bardelta - 2 + 2*exp(-barDelta) + 2*exp(-bardelta) ...
    - exp(-(barDelta-bardelta)) - exp(-(barDelta+bardelta)) );

ds = 1/4 * ( 2*dbardelta + 2*exp(-barDelta).*(-dbarDelta)...
    + 2*exp(-bardelta).*(-dbardelta) ...
    - exp(-(barDelta-bardelta)).*(-(dbarDelta-dbardelta))...
    - exp(-(barDelta+bardelta)).*(-(dbarDelta+dbardelta)) );

D = s*w^2*td^2./n.^2./(t-n/3);
dD1 = 2*w*td*s./n.^2./(t-n/3);
dD2 = (2*w^2*td*s + 2^2*td^2*ds)./n.^2./(t-n/3);
dD = [dD1, dD2];

end
