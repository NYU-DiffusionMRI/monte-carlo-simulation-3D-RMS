function D = vg(delta, Delta, r, D0)

        td=r^2/D0;
        bardelta=delta/td; dbardelta = -2*delta*D0 / r^3;
        barDelta=Delta/td; dbarDelta = -2*Delta*D0 / r^3;


        N=15; 
        b = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644 24.3113   27.4571   30.6019 ...
             33.7462   36.8900   40.0334   43.1766   46.3196 49.4624   52.6050   55.7476   58.8900   62.0323];
   
   
        s = 0; % ds = 0;
        for k=1:N
           s = s + (2/(b(k)^6*(b(k)^2-1)))*(-2 + 2*b(k)^2*bardelta + ...
                          2*(exp(-b(k)^2*bardelta)+exp(-b(k)^2*barDelta)) - ...
                          exp(-b(k)^2*(bardelta+barDelta)) - exp(-b(k)^2*(barDelta-bardelta))); 
%            ds = ds +  (2/(b(k)^6*(b(k)^2-1)))*( ...
%                             2*b(k)^2*dbardelta + ...
%                             2*exp(-b(k)^2*bardelta) - b(k)^2*dbardelta + ...
%                             2*exp(-b(k)^2*barDelta) - b(k)^2*dbarDelta - ...
%                             exp(-b(k)^2*(bardelta+barDelta)) - b(k)^2*(dbardelta + dbarDelta) - ... 
%                             exp(-b(k)^2*(barDelta-bardelta)) - b(k)^2*(dbarDelta - dbardelta));           
        end
        s = s.*D0.*td^3;
        D = s./(Delta-delta/3)./delta.^2;
%         ds = ds.*D0.*q.^2.*td^3 + 6*s.*q.^2.*r^5 / D0^2;
end