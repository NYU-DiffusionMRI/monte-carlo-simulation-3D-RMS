function D = un(t,n,w,td)

bardelta = n/td; % dbardelta = -n / td^3;
barDelta = t/td; % dbarDelta = -t / td^3;

s = 1/4 * ( 2*bardelta - 2 + 2*exp(-barDelta) + 2*exp(-bardelta) ...
    - exp(-(barDelta-bardelta)) - exp(-(barDelta+bardelta)) );

% ds = 1/4 * ( 2*dbardelta + 2*exp(-barDelta).*(-dbarDelta)...
%     + 2*exp(-bardelta).*(-dbardelta) ...
%     - exp(-(barDelta-bardelta)).*(-(dbarDelta-dbardelta))...
%     - exp(-(barDelta+bardelta)).*(-(dbarDelta+dbardelta)) );

D = s*w^2*td^2./n.^2./(t-n/3);
% dD1 = 2*w*td*s./n.^2./(t-n/3);
% dD2 = (2*w^2*td*s + 2^2*td^2*ds)./n.^2./(t-n/3);
% dD = [dD1, dD2];

end