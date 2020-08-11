% Part 2: Why elastic collision is the most reliable particle-membrane 
%    interaction?
% Demo 1: Equal-step-length random leap (Lee, et al., J Neurosci Methods
%    2020; Xing et al., MRM 2013)
%
% The purpose of this demontration includes:
% 1. Perform simple Monte Carlo simulations of diffusion between 1d, 2d, 
%    and 3d impermeable parallel planes, and show the bias in the 
%    diffusivity transverse and parallel to the membranes due to
%    inhomogeneous particle density caused by the ERL.
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
time_max = 10;                      % Maximal diffusion time (ms)
diffusivity = 2;                    % Intrinsic diffusivity (micron2/ms)
time_step = [0.02 0.04 0.06 0.08 0.1].^2/2/D0;      % Time for each step (ms)
membrane_distance = 1;              % Distance bewteen parallel planes (micron)

ia = interaction();
% Diffusion time, diffusivity transverse to membranes
t1d = zeros(1e3,numel(time_step)); D1d = zeros(1e3,numel(time_step));
t2d = zeros(1e3,numel(time_step)); D2d = zeros(1e3,numel(time_step));
t3d = zeros(1e3,numel(time_step)); D3d = zeros(1e3,numel(time_step));
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
    [t1d(:,i), D1d(:,i)] = ia.ERLdiffusivity1d(X);
    [t2d(:,i), D2d(:,i)] = ia.ERLdiffusivity2d(X);
    [t3d(:,i), D3d(:,i), AD3d(:,i)] = ia.ERLdiffusivity3d(X);
end
toc;
save(fullfile(root,'ERL_diffusivity.mat'),...
    't1d','t2d','t3d','D1d','D2d','D3d','AD3d',...
    'time_max','diffusivity','time_step','membrane_distance');

%% Plot diffusivity transverse to membranes
load(fullfile(root,'ERL_diffusivity.mat'));
D0 = diffusivity;
a = membrane_distance;

figure('unit','inch','position',[0 0 10 15]);

% Simulated diffusivity in 1d
subplot(321); hold on;
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
clear h
for i = 1:numel(time_step)
    h(i) = plot(1./t1d(:,i),D1d(:,i)/D0,'-','linewidth',1,'color',cmap(i,:));
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 1]); ylim([0 0.05]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','northwest');
pbaspect([1 1 1]);
box on; grid on;
title('1{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',20);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',20);

% Bias in distance between membranes in 1d
subplot(322)
hold on;
for i = 1:numel(time_step)
    D = D1d(:,i);
    t = t1d(:,i);
    [~,It] = min(abs(t-8));
    tlist = It:numel(t);
    aeff = sqrt((1/12./t(tlist))\D(tlist));
    plot(sqrt(2*D0*time_step(i))/4*2,a-aeff,'o','markersize',10,'linewidth',1); 
end
xlim([0 0.05]); ylim([0 0.05]);
set(gca,'fontsize',12);
hr = refline(1,0); set(hr,'color','k');
pbaspect([1 1 1]);
box on; grid on;
title('1{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$2\times\frac{1}{4}\delta s$ ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('$a-a^\prime$ ($\mu$m)','interpreter','latex','fontsize',20);

% Simulated diffusivity in 2d
subplot(323); hold on;
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
clear h
for i = 1:numel(time_step)
    h(i) = plot(1./t2d(:,i),D2d(:,i)/D0,'-','linewidth',1,'color',cmap(i,:));
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 1]); ylim([0 0.05]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','northwest');
pbaspect([1 1 1]);
box on; grid on;
title('2{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',20);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',20);

% Bias in distance between membranes in 2d
subplot(324)
hold on;
for i = 1:numel(time_step)
    D = D2d(:,i);
    t = t2d(:,i);
    [~,It] = min(abs(t-8));
    tlist = It:numel(t);
    aeff = sqrt((1/12./t(tlist))\D(tlist));
    plot(sqrt(4*D0*time_step(i))/2/pi*2,a-aeff,'o','markersize',10,'linewidth',1); 
end
xlim([0 0.05]); ylim([0 0.05]);
set(gca,'fontsize',12);
hr = refline(1,0); set(hr,'color','k');
pbaspect([1 1 1]);
box on; grid on;
title('2{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$2\times \frac{1}{2\pi}\,\delta s$ ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('$a-a^\prime$ ($\mu$m)','interpreter','latex','fontsize',20);

% Simulated diffusivity in 3d
subplot(325); hold on;
cmap = colormap('lines');
lgtxt = cell(numel(time_step),1);
clear h
for i = 1:numel(time_step)
    h(i) = plot(1./t3d(:,i),D3d(:,i)/D0,'-','linewidth',1,'color',cmap(i,:));
    lgtxt{i} = sprintf('$\\delta t = %.1f$ $\\mu$s',time_step(i)*1e3);
end
xlim([0 1]); ylim([0 0.05]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','northwest');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',20);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',20);

% Bias in distance between membranes in 3d
subplot(326)
hold on;
for i = 1:numel(time_step)
    D = D3d(:,i);
    t = t3d(:,i);
    [~,It] = min(abs(t-8));
    tlist = It:numel(t);
    aeff = sqrt((1/12./t(tlist))\D(tlist));
    plot(sqrt(6*D0*time_step(i))/8*2,a-aeff,'o','markersize',10,'linewidth',1); 
end
xlim([0 0.05]); ylim([0 0.05]);
set(gca,'fontsize',12);
hr = refline(1,0); set(hr,'color','k');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$2\times\frac{1}{8}\delta s$ ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('$a-a^\prime$ ($\mu$m)','interpreter','latex','fontsize',20);

%% Plot diffusivity parallel to membranes
load(fullfile(root,'ERL_diffusivity.mat'));
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
xlim([0 10]); ylim([1 1.05]);
set(gca,'fontsize',12);
legend(h,lgtxt,'interpreter','latex','fontsize',16,...
    'location','northeast');
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
    [~,It] = min(abs(t-5));
    tlist = It:numel(t);
    zeta(i) = mean(D(tlist))/D0 - 1;
    dx = sqrt(6*D0*time_step(i));
    SV = 2/a;
    plot(1/16*(SV*dx)./(1-1/4*SV*dx),zeta(i),'o','markersize',10,'linewidth',1); 
end
xlim([0 0.025]); ylim([0 0.025]);
set(gca,'fontsize',12);
hr = refline(1,0); set(hr,'color','k');
pbaspect([1 1 1]);
box on; grid on;
title('3{\itd} parallel planes','fontsize',20,'fontweight','normal');
xlabel('$\left(\frac{S}{V}\right)\cdot\delta s/16$','interpreter','latex','fontsize',20);
ylabel('$D_\parallel/D_0-1$','interpreter','latex','fontsize',20);

