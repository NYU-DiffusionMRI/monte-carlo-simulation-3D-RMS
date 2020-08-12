% Part 2: Why elastic collision is the most reliable particle-membrane 
%    interaction?
% Demo 2: Equal-step-length random leap (Lee, et al., J Neurosci Methods
%    2020; Ford and Hackney 1997)
%
% The purpose of this demontration includes:
% 2. Perform simple Monte Carlo simulations of diffusion between 1d, 2d, 
%    and 3d permeable parallel planes to demonstrate the compatibility of 
%    rejection sampling with the simulation of water exchange.
%
% Reference: 
% Lee, et al., Journal of Neuroscience Methods 2020 (in preparation)
% Ford and Hackney, Magnetic Resonance in Medicine 1997 (doi:10.1002/mrm.1910370315)
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

% Setup directory to this code
root = pwd;

% Load simulation tool
addpath(fullfile(root,'analysis'));

%% Simulation of particle density around permeable membranes

% Diffusion between permeable parallel planes in 1d, 2d and 3d
% Permeable parallel planes are placed at x = +a/2 and -a/2
% Particles are initialized at x = 0
% Particle density is recorded between x = 0 and a
time_max = 0.2;                     % Maximal diffusion time (ms)
diffusivity = 0.5;                  % Intrinsic diffusivity (micron2/ms)
time_step = 0.02.^2/2/diffusivity;  % Time for each step (ms)
bin_num = 1000;                     % # bin to calculate particle density
membrane_distance = 1;              % Distance bewteen parallel planes (micron)
permeability = 0.1;                 % Permeability (micron/ms)

ia = interaction();
tic;
X.particle_num = 1e7;               % # particle
X.time_max = time_max;              % Maximal diffusion time (ms)
X.diffusivity = diffusivity;        % Intrinsic diffusivity (micron2/ms)
X.time_step = time_step;            % Time for each step (ms)
X.bin_num = bin_num;                % # bin to calculate particle density
% Distance bewteen parallel planes (micron)
X.membrane_distance = membrane_distance;
X.permeability = permeability;      % Permeability (micron/ms)
[t1d,density1d] = ia.RSpermeability1d(X);
[t2d,density2d] = ia.RSpermeability2d(X);
[t3d,density3d] = ia.RSpermeability3d(X);
toc;
save(fullfile(root,'RS_permeability.mat'),...
    't1d','t2d','t3d','density1d','density2d','density3d',...
    'time_max','diffusivity','time_step','membrane_distance','permeability');

%% Plot permeability based on the simulated particle density
% Applying > 1e7 particles is strongly suggested for a reliable result
load(fullfile(root,'RS_permeability.mat'));

figure('unit','inch','position',[0 0 15 5]);

% Permeability in 1d
subplot(131); hold on;
fit_num = 3;
ia = interaction();
% The particle density is discretized in 1d. Resize the matrix using box
% can solve the problem.
density_resize = imresize(density1d,size(density1d).*[1 50/1000],'box');
kappa1d = ia.permeability(density_resize,membrane_distance,diffusivity,fit_num);
plot(t1d,kappa1d,'-','linewidth',0.75);
xlim([0 time_max]); ylim([0 0.5]);
hr = refline(0,permeability); set(hr,'color','r','linewidth',1);
legend(hr,'Ground truth','fontsize',16);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('1{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$\kappa$ ($\mu$m/ms)','interpreter','latex','fontsize',20);

% Permeability in 2d
subplot(132); hold on;
fit_num = 50;
ia = interaction();
kappa2d = ia.permeability(density2d,membrane_distance,diffusivity,fit_num);
plot(t2d,kappa2d,'-','linewidth',0.75);
xlim([0 time_max]); ylim([0 0.5]);
hr = refline(0,permeability); set(hr,'color','r','linewidth',1);
legend(hr,'Ground truth','fontsize',16);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('2{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$\kappa$ ($\mu$m/ms)','interpreter','latex','fontsize',20);

% Permeability in 3d
subplot(133); hold on;
fit_num = 50;
ia = interaction();
kappa3d = ia.permeability(density3d,membrane_distance,diffusivity,fit_num);
plot(t3d,kappa3d,'-','linewidth',0.75);
xlim([0 time_max]); ylim([0 0.5]);
hr = refline(0,permeability); set(hr,'color','r','linewidth',1);
legend(hr,'Ground truth','fontsize',16);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$\kappa$ ($\mu$m/ms)','interpreter','latex','fontsize',20);


