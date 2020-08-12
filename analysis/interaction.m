classdef interaction < handle
    
    properties (GetAccess = public, SetAccess = public)
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
    end
    
    methods (Access = public)
        function this = interaction()
            
        end
        
        function density = ERLdensity1d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(2*D0*dt);         % Step size (micron)
            
            xi = a*rand(Np,1);          % Initial position (micron)
            density = zeros(Nbin,1);
            for i = 1:Np
                % Step in random direction
                dxi = (2*(rand(Nt,1)>0.5)-1)*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        xk = xj + (2*(rand>0.5)-1)*dx;
                    end
                    xj = xk;
                    if j == Nt          % Calculate particle density
                        k = ceil(xj/(a/Nbin));
                        density(k) = density(k) + 1;
                    end
                end
            end
        end
        
        function density = ERLdensity2d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(4*D0*dt);         % Step size (micron)
            
            xi = a*rand(Np,1);          % Initial position (micron)
            density = zeros(Nbin,1);
            for i = 1:Np
                % Step in random direction
                dxi = cos(2*pi*rand(Nt,1))*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        xk = xj + cos(2*pi*rand)*dx;
                    end
                    xj = xk;
                    if j == Nt          % Calculate particle density
                        k = ceil(xj/(a/Nbin));
                        density(k) = density(k) + 1;
                    end
                end
            end
        end
        
        function density = ERLdensity3d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(6*D0*dt);         % Step size (micron)
            
            xi = a*rand(Np,1);          % Initial position (micron)
            density = zeros(Nbin,1);
            for i = 1:Np
                % Step in random direction
                dxi = (2*rand(Nt,1)-1)*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        xk = xj + (2*rand-1)*dx;
                    end
                    xj = xk;
                    if j == Nt          % Calculate particle density
                        k = ceil(xj/(a/Nbin));
                        density(k) = density(k) + 1;
                    end
                end
            end
        end
        
        function h = theory_ERLdensity1d(this,diffusivity,time_step,bin_num,membrane_distance)
            % Impermeable parallel planes are placed at x = 0 and a
            D0 = diffusivity;
            dt = time_step;
            Nbin = bin_num;
            a = membrane_distance;
            
            xx = linspace(0,a,Nbin+1)/a;
            xx = xx(1:end-1)/2 + xx(2:end)/2;
            yy = ones(Nbin,1);
            dx = sqrt(2*D0*dt);
            xs = floor(dx/a*Nbin);
            yy(1:xs) = 0.5;
            yy(end-xs+1:end) = 0.5;
            yy = yy/mean(yy);
            h = plot(xx,yy,'-');
        end
        
        function h = theory_ERLdensity2d(this,diffusivity,time_step,bin_num,membrane_distance)
            % Impermeable parallel planes are placed at x = 0 and a
            D0 = diffusivity;
            dt = time_step;
            Nbin = bin_num;
            a = membrane_distance;
            
            xx = linspace(0,a,Nbin+1)/a;
            xx = xx(1:end-1)/2 + xx(2:end)/2;
            yy = ones(Nbin,1);
            dx = sqrt(4*D0*dt);
            xs = floor(dx/a*Nbin);
            yy(1:xs) = 1-acos(xx(1:xs)*a/dx)/pi;
            yy(end-xs+1:end) = 1-acos(xx(xs:-1:1)*a/dx)/pi;
            yy = yy/mean(yy);
            h = plot(xx,yy,'-');
        end
        
        function h = theory_ERLdensity3d(this,diffusivity,time_step,bin_num,membrane_distance)
            % Impermeable parallel planes are placed at x = 0 and a
            D0 = diffusivity;
            dt = time_step;
            Nbin = bin_num;
            a = membrane_distance;
            
            xx = linspace(0,a,Nbin+1)/a;
            xx = xx(1:end-1)/2 + xx(2:end)/2;
            yy = ones(Nbin,1);
            dx = sqrt(6*D0*dt);
            xs = floor(dx/a*Nbin);
            yy(1:xs) = 1/2+(xx(1:xs)*a/dx)/2;
            yy(end-xs+1:end) = 1/2+(xx(xs:-1:1)*a/dx)/2;
            yy = yy/mean(yy);
            h = plot(xx,yy,'-');
        end
        
        function h = plotdensity(this,density,membrane_distance,scale)
            Nbin = size(density,2);
            a = membrane_distance;
            
            xx = linspace(0,a,Nbin+1)/a;
            xx = xx(1:end-1)/2 + xx(2:end)/2;
            xx = xx(1:end/2);
            yy = density(:).';
            yy = yy(1:end/2)/2 + yy(end:-1:end/2+1)/2;
            yy = yy/(a/Nbin)/sum(density);
            xx = imresize(xx,scale,'box');
            yy = imresize(yy,scale,'box');
            h = plot(xx,yy,'.-');
        end
        
        function [t,D] = ERLdiffusivity1d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Ni = 1e3;                   % # recorded diffusion time
            Nt = ceil(Tmax/dt/Ni)*Ni;   % # step
            dx = sqrt(2*D0*dt);         % Step size (micron)
            Tstep = Nt/Ni;              % # step between recorded time
            
            x2 = zeros(Ni,1);           % Second order cumulant of displacement
            xi = a*rand(Np,1);          % Initial position (micron)
            for i = 1:Np
                dxi = (2*(rand(Nt,1)>0.5)-1)*dx;
                xj = xi(i);
                k = 1;
                for j = 1:Nt
                    xk = xj + dxi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        xk = xj + (2*(rand>0.5)-1)*dx;
                    end
                    xj = xk;
                    if j == (Tstep*k)
                        x2(k) = x2(k) + (xj-xi(i))*(xj-xi(i));
                        k = k+1;
                    end
                end
            end
            t = (Tstep:Tstep:Nt)*dt; t = t(:);
            D = x2/2./t/Np;
        end
        
        function [t,D] = ERLdiffusivity2d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Ni = 1e3;                   % # recorded diffusion time
            Nt = ceil(Tmax/dt/Ni)*Ni;   % # step
            dx = sqrt(4*D0*dt);         % Step size (micron)
            Tstep = Nt/Ni;              % # step between recorded time
            
            x2 = zeros(Ni,1);           % Second order cumulant of displacement
            xi = a*rand(Np,1);          % Initial position (micron)
            for i = 1:Np
                dxi = cos(2*pi*rand(Nt,1))*dx;
                xj = xi(i);
                k = 1;
                for j = 1:Nt
                    xk = xj + dxi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        xk = xj + cos(2*pi*rand)*dx;
                    end
                    xj = xk;
                    if j == (Tstep*k)
                        x2(k) = x2(k) + (xj-xi(i))*(xj-xi(i));
                        k = k+1;
                    end
                end
            end
            t = (Tstep:Tstep:Nt)*dt; t = t(:);
            D = x2/2./t/Np;
        end
        
        function [t,RD,AD] = ERLdiffusivity3d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Ni = 1e3;                   % # recorded diffusion time
            Nt = ceil(Tmax/dt/Ni)*Ni;   % # step
            dx = sqrt(6*D0*dt);         % Step size (micron)
            Tstep = Nt/Ni;              % # step between recorded time
            
            x2 = zeros(Ni,1);           % Second order cumulant of displacement transverse to membrane
            y2 = zeros(Ni,1);           % Second order cumulant of displacement parallel to membrane
            xi = a*rand(Np,1);          % Initial position (micron)
            for i = 1:Np
                cos_theta = (2*rand(Nt,1)-1);
                dxi = cos_theta*dx;
                dyi = sqrt(1-cos_theta.^2).*sin(2*pi*rand(Nt,1))*dx;
                xj = xi(i);
                yj = 0;
                k = 1;
                for j = 1:Nt
                    xk = xj + dxi(j);
                    yk = yj + dyi(j);
                    while xk<0 || xk>a  % Equal-step-length random leap
                        cos_thetai = (2*rand-1);
                        xk = xj + cos_thetai*dx;
                        yk = yj + sqrt(1-cos_thetai^2)*sin(2*pi*rand)*dx;
                    end
                    xj = xk;
                    yj = yk;
                    if j == (Tstep*k)
                        x2(k) = x2(k) + (xj-xi(i))*(xj-xi(i));
                        y2(k) = y2(k) + yj*yj;
                        k = k+1;
                    end
                end
            end
            t = (Tstep:Tstep:Nt)*dt; t = t(:);
            RD = x2/2./t/Np;
            AD = y2/2./t/Np;
        end
        
        function [t,RD,AD] = RSdiffusivity3d(this,X)
            % Impermeable parallel planes are placed at x = 0 and a
            % Particles are homogeneously initialized between x = 0 and a
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            
            Ni = 1e3;                   % # recorded diffusion time
            Nt = ceil(Tmax/dt/Ni)*Ni;   % # step
            dx = sqrt(6*D0*dt);         % Step size (micron)
            Tstep = Nt/Ni;              % # step between recorded time
            
            x2 = zeros(Ni,1);           % Second order cumulant of displacement transverse to membrane
            y2 = zeros(Ni,1);           % Second order cumulant of displacement parallel to membrane
            xi = a*rand(Np,1);          % Initial position (micron)
            for i = 1:Np
                cos_theta = (2*rand(Nt,1)-1);
                dxi = cos_theta*dx;
                dyi = sqrt(1-cos_theta.^2).*sin(2*pi*rand(Nt,1))*dx;
                xj = xi(i);
                yj = 0;
                k = 1;
                for j = 1:Nt
                    xk = xj + dxi(j);
                    yk = yj + dyi(j);
                    if xk<0 || xk>a  % Rejection sampling
                        xk = xj;
                        yk = yj;
                    end
                    xj = xk;
                    yj = yk;
                    if j == (Tstep*k)
                        x2(k) = x2(k) + (xj-xi(i))*(xj-xi(i));
                        y2(k) = y2(k) + yj*yj;
                        k = k+1;
                    end
                end
            end
            t = (Tstep:Tstep:Nt)*dt; t = t(:);
            RD = x2/2./t/Np;
            AD = y2/2./t/Np;
        end
        
        function [t,density] = RSpermeability1d(this,X)
            % Permeable parallel planes are placed at x = +a/2 and -a/2
            % Particles are initialized at x = 0
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            kappa = X.permeability;     % Permeability (micron/ms)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(2*D0*dt);         % Step size (micron)
            Pm = kappa*dx/D0;           % Permeation probability
            Pm = Pm/(1+Pm);             % Correction for the permeation probability
            
            xi = zeros(Np,1);           % Initial position (micron)
            density = zeros(Nt,Nbin);
            for i = 1:Np
                % Step in random direction
                dxi = (2*(rand(Nt,1)>0.5)-1)*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    % Rejection sampling
                    if ( (xj+a/2)*(xk+a/2)<0 || (xj-a/2)*(xk-a/2)<0 ) && rand<(1-Pm)
                        xk = xj;
                    end
                    xj = xk;
                    % Calculate particle density between x = 0 and a
                    if abs(xj)>0 && abs(xj)<a
                        k = ceil(abs(xj)/(a/Nbin));
                        density(j,k) = density(j,k) + 1;
                    end
                end
            end
            t = (1:Nt)*dt; t = t(:);
        end
        
        function [t,density] = RSpermeability2d(this,X)
            % Permeable parallel planes are placed at x = +a/2 and -a/2
            % Particles are initialized at x = 0
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            kappa = X.permeability;     % Permeability (micron/ms)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(4*D0*dt);         % Step size (micron)
            Pm = pi/4*kappa*dx/D0;      % Permeation probability
            Pm = Pm/(1+Pm);             % Correction for the permeation probability
            
            xi = zeros(Np,1);           % Initial position (micron)
            density = zeros(Nt,Nbin);
            for i = 1:Np
                % Step in random direction
                dxi = cos(2*pi*rand(Nt,1))*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    % Rejection sampling
                    if ( (xj+a/2)*(xk+a/2)<0 || (xj-a/2)*(xk-a/2)<0 ) && rand<(1-Pm)
                        xk = xj;
                    end
                    xj = xk;
                    % Calculate particle density between x = 0 and a
                    if abs(xj)>0 && abs(xj)<a
                        k = ceil(abs(xj)/(a/Nbin));
                        density(j,k) = density(j,k) + 1;
                    end
                end
            end
            t = (1:Nt)*dt; t = t(:);
        end
        
        function [t,density] = RSpermeability3d(this,X)
            % Permeable parallel planes are placed at x = +a/2 and -a/2
            % Particles are initialized at x = 0
            Np = X.particle_num;        % # particle
            Tmax = X.time_max;          % Maximal diffusion time (ms)
            D0 = X.diffusivity;         % Intrinsic diffusivity (micron2/ms)
            dt = X.time_step;           % Time for each step (ms)
            Nbin = X.bin_num;           % # bin to calculate particle density
            a = X.membrane_distance;    % Distance bewteen parallel planes (micron)
            kappa = X.permeability;     % Permeability (micron/ms)
            
            Nt = ceil(Tmax/dt);         % # step
            dx = sqrt(6*D0*dt);         % Step size (micron)
            Pm = 2/3*kappa*dx/D0;       % Permeation probability
            Pm = Pm/(1+Pm);             % Correction for the permeation probability
            
            xi = zeros(Np,1);           % Initial position (micron)
            density = zeros(Nt,Nbin);
            for i = 1:Np
                % Step in random direction
                dxi = (2*rand(Nt,1)-1)*dx;
                xj = xi(i);             % Initial position
                for j = 1:Nt
                    xk = xj + dxi(j);
                    % Rejection sampling
                    if ( (xj+a/2)*(xk+a/2)<0 || (xj-a/2)*(xk-a/2)<0 ) && rand<(1-Pm)
                        xk = xj;
                    end
                    xj = xk;
                    % Calculate particle density between x = 0 and a
                    if abs(xj)>0 && abs(xj)<a
                        k = ceil(abs(xj)/(a/Nbin));
                        density(j,k) = density(j,k) + 1;
                    end
                end
            end
            t = (1:Nt)*dt; t = t(:);
        end
        
        function kappa = permeability(this,density,membrane_distance,diffusivity,fit_num)
            % Permeable membrane is placed at the center
            % Particle is initialized over the left side of the membrane
            [Nt,Nbin] = size(density);
            a = membrane_distance;
            D0 = diffusivity;
            
            kappa = zeros(Nt,1);
            for i = 1:Nt
                xx = linspace(0,a,Nbin+1)/a;
                xx = xx(1:end-1)/2 + xx(2:end)/2;
                yy = density(i,:).';
                yy = yy/(a/Nbin)/sum(yy);
                xx = xx(:); yy = yy(:);
                nfit = fit_num;
                [~,Im] = min(abs(xx-0.5));
                llist = Im-nfit:Im-1;
                rlist = Im+1:Im+nfit;
                ltmp = [ones(nfit,1), xx(llist)-0.5]\yy(llist);
                rtmp = [ones(nfit,1), xx(rlist)-0.5]\yy(rlist);
                kappa(i) = D0*rtmp(2)/(rtmp(1)-ltmp(1));
            end
        end
        
        
    end
    
        
end