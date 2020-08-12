% Part 2: Why elastic collision is the most reliable particle-membrane 
%    interaction?
% Demo 2: Equal-step-length random leap (Lee, et al., J Neurosci Methods
%    2020; Ford and Hackney 1997)
%
% The purpose of this demontration includes:
% 1. Perform simple Monte Carlo simulations of diffusion between 1d, 2d, 
%    and 3d impermeable parallel planes to show the bias in the diffusivity
%    parallel to membranes caused by rejection sampling.
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

%% Simulation of diffusivity parallel to membranes

% Diffusion between impermeable parallel planes in 1d, 2d and 3d
% Impermeable parallel planes are placed at x = 0 and a
% Particles are initialized homogeneously between x = 0 and a
time_max = 10;                      % Maximal diffusion time (ms)
diffusivity = 2;                    % Intrinsic diffusivity (micron2/ms)
time_step = [0.02 0.04 0.06 0.08 0.1].^2/2/diffusivity;      % Time for each step (ms)
membrane_distance = 1;              % Distance bewteen parallel planes (micron)

ia = interaction();
% Diffusion time
t3d = zeros(1e3,numel(time_step));
% Diffusivity parallel to membranes in 3d
AD3d = zeros(1e3,numel(time_step));
tic;
for i = 1:numel(time_step)
    X.particle_num = 1e5;           % # particle
    X.time_max = time_max;          % Maximal diffusion time (ms)
    X.diffusivity = diffusivity;    % Intrinsic diffusivity (micron2/ms)
    X.time_step = time_step(i);     % Time for each step (ms)
    % Distance bewteen parallel planes (micron)
    X.membrane_distance = membrane_distance;
    [t3d(:,i), ~, AD3d(:,i)] = ia.RSdiffusivity3d(X);
end
toc;
save(fullfile(root,'RS_diffusivity.mat'),...
    't3d','AD3d','time_max','diffusivity','time_step','membrane_distance');

%% Plot diffusivity parallel to membranes
% Applying > 1e6 particles is strongly suggested for a reliable result
load(fullfile(root,'RS_diffusivity.mat'));
D0 = diffusivity;
a = membrane_distance;

figure('unit','inch','position',[0 0 10 5]);

% Simulated diffusivity parallel to membranes in 3d
subplot(121); hold on;
lgtxt = cell(numel(time_step),1);
clear h
for i = 1:numel(time_step)
    h(i) = plot(t3d(:,i),AD3d(:,i)/D0,'-','linewidth',1);
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 10]); ylim([0.9 1]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','south');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$D_\parallel(t)/D_0$','interpreter','latex','fontsize',20);

% Bias in diffusivity parallel to membranes in 3d
subplot(122)
hold on;
zeta = zeros(numel(time_step),1);
for i = 1:numel(time_step)
    D = AD3d(:,i);
    t = t3d(:,i);
    [~,It] = min(abs(t-8));
    tlist = It:numel(t);
    zeta(i) = mean(D(tlist))/D0 - 1;
    dx = sqrt(6*D0*time_step(i));
    SV = 2/a;
    plot(3/16*(SV*dx),zeta(i),'o','markersize',10,'linewidth',1); 
end
xlim([0 0.08]); ylim([-0.08 0]);
set(gca,'fontsize',12,'ytick',-0.08:0.02:0);
hr = refline(-1,0); set(hr,'color','k');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$\left(\frac{S}{V}\right)\cdot\frac{3}{16}\delta s$','interpreter','latex','fontsize',20);
ylabel('$D_\parallel/D_0-1$','interpreter','latex','fontsize',20);
