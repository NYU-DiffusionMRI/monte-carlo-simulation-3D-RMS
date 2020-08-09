function X = sticksizefit_app(S,b)

S = double(S(:)); b = double(b(:));

Nt = 1e1;
Xt = zeros(2,Nt);
rng(0);
for i = 1:Nt
    x0 = [sqrt(2); 0.05+0.1*rand];
    options = optimoptions(@lsqnonlin,'Display','off','Algorithm','Levenberg-Marquardt');
    Xt(:,i) = lsqnonlin(@(x)costfunction(S,b,x),x0,[],[],options);
end
% idx = dbscan(Xt.',0.2,3);
idx = kmeans(Xt.',3);
Xt = Xt(:,idx==mode(idx));
X = zeros(2,1);
X(1) = median(Xt(1,:));
[~,I] = min(abs( Xt(1,:)-median(Xt(1,:)) ));
X(1) = X(1)^2;
X(2) = Xt(2,I);

end

function [J, dJ] = costfunction(S,b,X)

nb = numel(b);
SX = stickmodel(b,X);
J = 1/2/nb*sum((SX-S).^2);

dS = zeros(nb,2);
dS(:,1) = -sqrt(pi/4./b).*exp(-b.*X(2))./X(1).^2;
dS(:,2) = sqrt(pi/4./b)./X(1).*exp(-b.*X(2)).*(-b);
% dS(:,1) = sqrt(pi/4./b).*exp(-b.*X(3)/scale);
% dS(:,2) = 0;
% dS(:,3) = X(1).*sqrt(pi/4./b).*exp(-b.*X(3)/scale).*(-b/scale);

dJ = zeros(2,1);
dJ(1) = 1/nb*sum((SX-S).*dS(:,1));
dJ(2) = 1/nb*sum((SX-S).*dS(:,2));

end

function S = stickmodel(b,X)

S = sqrt(pi/4./b)./X(1).*exp(-b.*X(2));

end
