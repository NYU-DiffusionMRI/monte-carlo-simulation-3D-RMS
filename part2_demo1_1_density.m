% Part 2: Why elastic collision is the most reliable particle-membrane 
%    interaction?
% Demo 1: Equal-step-length random leap (Lee, et al., J Neurosci Methods
%    2020; Xing et al., MRM 2013)
%
% The purpose of this demontration includes:
% 1. Perform simple Monte Carlo simulations of diffusion between 1d, 2d, 
%    and 3d impermeable parallel planes, and show the inhomogeneous 
%    particle density around membranes caused by the ERL.
%
% Reference: 
% Lee, et al., Journal of Neuroscience Methods 2020 (in preparation)
% Xing, et al., Magnetic Resonance in Medicine 2013 (doi:10.1002/mrm.24551)
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

% Setup directory to this code
root = pwd;

% Load simulation tool
addpath(fullfile(root,'analysis'));

%% Simulation of particle density around impermeable membranes

% Diffusion between impermeable parallel planes in 1d, 2d and 3d
% Impermeable parallel planes are placed at x = 0 and a
% Particles are initialized homogeneously between x = 0 and a
% Particle density is recorded between x = 0 and a
time_max = 10;                      % Maximal diffusion time (ms)
diffusivity = 2;                    % Intrinsic diffusivity (micron2/ms)
time_step = [0.02 0.04 0.06 0.08 0.1].^2/2/diffusivity;      % Time for each step (ms)
bin_num = 1000;                     % # bin to calculate particle density
membrane_distance = 1;              % Distance bewteen parallel planes (micron)

ia = interaction();
density1d = zeros(numel(time_step),bin_num);
density2d = zeros(numel(time_step),bin_num);
density3d = zeros(numel(time_step),bin_num);
tic;
for i = 1:numel(time_step)
    X.particle_num = 1e5;           % # particle
    X.time_max = time_max;          % Maximal diffusion time (ms)
    X.diffusivity = diffusivity;    % Intrinsic diffusivity (micron2/ms)
    X.time_step = time_step(i);     % Time for each step (ms)
    X.bin_num = bin_num;            % # bin to calculate particle density
    % Distance bewteen parallel planes (micron)
    X.membrane_distance = membrane_distance;
    density1d(i,:) = ia.ERLdensity1d(X);
    density2d(i,:) = ia.ERLdensity2d(X);
    density3d(i,:) = ia.ERLdensity3d(X);
end
toc;
save(fullfile(root,'ERL_density.mat'),...
    'density1d','density2d','density3d',...
    'time_max','diffusivity','time_step','membrane_distance');

%% Plot particle density around impermeable membranes
% Applying > 1e6 particles is strongly suggested for a reliable result
load(fullfile(root,'ERL_density.mat'));

figure('unit','inch','position',[0 0 10 15]);

% Theory in 1d
subplot(321); hold on;
bin_num = size(density1d,2);
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.theory_ERLdensity1d(diffusivity,time_step(i),bin_num,membrane_distance);
    set(h(i),'linewidth',1,'color',cmap(i,:));
    hr = xline(sqrt(2*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','southeast');
pbaspect([1 1 1]);
box on; grid on;
title('1{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);

% Simulation in 1d
subplot(322);
hold on;
scale = 0.25;                       % Scale to resize the density array
cmap = colormap('lines');
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.plotdensity(density1d(i,:),membrane_distance,scale);
    set(h(i),'linewidth',1,'color',cmap(i,:),'markersize',10);
    hr = xline(sqrt(2*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('1{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);

% Theory in 2d
subplot(323); hold on;
bin_num = size(density2d,2);
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.theory_ERLdensity2d(diffusivity,time_step(i),bin_num,membrane_distance);
    set(h(i),'linewidth',1,'color',cmap(i,:));
    hr = xline(sqrt(4*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','southeast');
pbaspect([1 1 1]);
box on; grid on;
title('2{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);

% Simulation in 2d
subplot(324);
hold on;
scale = 0.25;                       % Scale to resize the density array
cmap = colormap('lines');
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.plotdensity(density2d(i,:),membrane_distance,scale);
    set(h(i),'linewidth',1,'color',cmap(i,:),'markersize',10);
    hr = xline(sqrt(4*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('2{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);

% Theory in 3d
subplot(325); hold on;
bin_num = size(density3d,2);
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.theory_ERLdensity3d(diffusivity,time_step(i),bin_num,membrane_distance);
    set(h(i),'linewidth',1,'color',cmap(i,:));
    hr = xline(sqrt(6*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','southeast');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);

% Simulation in 3d
subplot(326);
hold on;
scale = 0.25;                       % Scale to resize the density array
cmap = colormap('lines');
ia = interaction();
clear h
for i = 1:numel(time_step)
    h(i) = ia.plotdensity(density3d(i,:),membrane_distance,scale);
    set(h(i),'linewidth',1,'color',cmap(i,:),'markersize',10);
    hr = xline(sqrt(6*diffusivity*time_step(i)));
    set(hr,'color',cmap(i,:),'linestyle','--','linewidth',1);
end
xlim([0 0.3]); ylim([0 1.2]);
set(gca,'fontsize',12);
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$h/a$','interpreter','latex','fontsize',20);
ylabel('$\tilde{\rho}(h)$','interpreter','latex','fontsize',20);
