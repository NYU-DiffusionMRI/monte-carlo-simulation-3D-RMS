% Part 1: Diffusivity time-dependence along realistic white matter axons
% Demo 4: Analysis (Lee, et al., Commun Biol 2020)
%
% The purpose of this demontration includes:
% 1. Calculate diffusivity and kurtosis time-dependence based on
%    displacement cumulants (Fig. 2b-h).
%
% Reference: 
% Lee, et al., Communications Biology 2020 (doi:10.1038/s42003-020-1050-x)
% Lee, et al., Brain Structure and Function 2019 
%    (doi:10.1007/s00429-019-01844-6)
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

% Setup directory to this code
root = pwd;

% Load RMS analysis tool
addpath(fullfile(root,'analysis'));

% Load SphereDistributionRand tool
addpath(genpath(fullfile(root,'tool')));

%% Read simulation results
target = cell(4,1);
target{1} = fullfile(root,'hpc_code','input','I_realistic_axon_mitochondria');
target{2} = fullfile(root,'hpc_code','input','II_realistic_axon');
target{3} = fullfile(root,'hpc_code','input','III_caliber_variation');
target{4} = fullfile(root,'hpc_code','input','IV_undulation');

rr = readRMS();
sim = struct([]);
for i = 1:numel(target)
    files = dir(fullfile(target{i},'setup*'));
    for j = 1:numel(files)
        sim(i,j).data = readRMS(files(j).name);
    end
end
save(fullfile(root,'hpc_code','result','example.mat'),'sim');

%% Analyze packing
fiber = struct([]);
for i = 1:numel(target)
    file = load(fullfile(target{i},'fiber.mat'));
    fiberi = file.fiber;
    for j = 1:numel(fiberi)
        ri = sqrt(sum(sum(fiberi(j).shape,1),2)/pi)*fiberi(j).voxel_size;
        fiber(i,j).CV_radius = std(ri)/mean(ri);
        fiber(i,j).volume = nnz(fiberi(j).shape)*fiberi(j).voxel_size^3;
        if i == 1
            fiber(i,j).mitochondiral_volume = ...
                nnz(fiberi(j).shape==2)*fiberi(j).voxel_size^3;
        end
    end
end

for i = 1:numel(target)
    volume = sum([fiber(i,:).volume]);
    for j = 1:numel(fiberi)
        fiber(i,j).volume_fraction = fiber(i,j).volume/volume;
    end
end

%% Calculate diffusivity and kurtosis based on displacement cumulants
load(fullfile(root,'hpc_code','result','example.mat'));

% Plot Fig. 2b-c: AD and AK in four geometries in Fig. 2a (no dispersion)
figure('unit','inch','position',[0 0 12 8]);
cmap = colormap('lines');
clear h
for i = 1:size(sim,1)
    Di = []; Kj = [];
    for j = 1:size(sim,2)
        simj = sim(i,j).data;
        n = [0 0 1];
        [Kj, Dj] = simj.akc_mom(n);
        Di = cat(2,Di,Dj);
        Ki = cat(2,Ki,Kj);
    end
    fi = [fiber(i,:).volume_fraction];
    D = sum(fi.*Di,2);
    K = 3*( sum(fi.*Di.^2,2) - D.^2 ) ./ D.^2 + sum(fi.*Di.^2.*Ki)./D.^2;
    t = simj.TD;
    
    subplot(121);
    hold on;
    plot(t,D,'.','MarkerSize',8,'Color',cmap(i,:));
    
    subplot(121);
    hold on;
    h(i) = plot(t,K,'.','MarkerSize',8,'Color',cmap(i,:));
end

subplot(121);
xlim([0 0.5]); ylim([0.6 0.95]);
set(gca,'fontsize',14);
xlabel('$1/\sqrt{t}$ (ms$^{-1/2}$)','Interpreter','latex','FontSize',24);
ylabel('$D(t)/D_0$','Interpreter','latex','FontSize',24);
pbaspect([1 1.5 1]); box on;

subplot(122);
xlim([0 0.5]); ylim([0 0.5]);
set(gca,'fontsize',14);
xlabel('$1/\sqrt{t}$ (ms$^{-1/2}$)','Interpreter','latex','FontSize',24);
ylabel('$K(t)$','Interpreter','latex','FontSize',24);
pbaspect([1 1.5 1]); box on;
legend(h,{'I: original','II: mitochondria removed',...
    'III. caliber variation only','IV: undulation only'},...
    'Interpreter','latex','FontSize',16,'Location','east');

% Plot Fig. 2d: Relation of relative diffusivity and coefficient of
% variation of radius
Da = sim(1,1).data.D0(2);       % Intrinsic diffusivity in cytoplasm
Dm = sim(1,1).data.D0(3);       % Intrinsic diffusivity in mitochondria
% Mitochondrial volume fraction
fm = [fiber(1,:).mitochondrial_volume]./[fiber(1,:).volume];
D0 = (1-fm)*Da + fm*Dm;         % Mean instrinsic diffusivity
zeta = D0./Dinf(1,:)-1;         % Relative diffusivity change
CV = [fiber(1,:).CV_radius];    % Coefficient of variation of radius

figure;
plot(CV.^2,zeta,'.','MarkerSize',10);
xlim([0 0.4]); ylim([0 3]);
set(gca,'fontsize',16);
xlabel('${\rm CV}^2(r)$','Interpreter','latex','FontSize',24);
ylabel('$\zeta$','Interpreter','latex','FontSize',24);
pbaspect([1 1 1]); box on; grid on;

% Plot Fig. e-f: AD and AK in geometry I with four dispersion angles
mu = [0 0 1];                   % Mean direction: z-axis
kappa = [Inf 15.4 4.7 1.65];    % polar angle: 0 ,15, 30, 45 degree
Ntheta = numel(kappa);          % # polar angle
Nphi = 10;                      % # azimuthal angle

simi = sim(1,:);
figure;
cmap = colormap('lines');
cmap = cmap([1 5:7],:);
clear h
for i = 1:Ntheta
    % Sampling directions based on the Watson distribution
    rng(0);
    if k == 1
        n = repmat([0 0 1],numel(simi)*Nphi,1);
    else
        n = randWatson(numel(simi)*Nphi, mu, kappa(i));
    end
    
    % Calculate displacement cumulants
    dx2i = []; dx4i = [];
    for j = 1:numel(simi)
        list = (j-1)*Nphi+1 : j*Nphi;
        ni = n(list,:);
        [Kj,Dj] = simi(j).data.akc_mom(ni);
        t = simi(j).data.TD;
        dx2j = Dj*2.*t;
        dx4j = (Kj+3).*dx2j.^2;
        dx2j = mean(dx2j,2);
        dx4j = mean(dx4j,2);
        dx2i = cat(2,dx2i,dx2j);
        dx4i = cat(2,dx4i,dx4j);
    end
    dx2 = sum(fi.*dx2i,2);
    dx4 = sum(fi.*dx4i,2);
    
    % Calculate diffusivity and kurtosis
    D = dx2/2./t;
    K = dx4./dx2.^2-3;
    
    subplot(121);
    hold on;
    plot(t,D,'.','MarkerSize',8,'Color',cmap(i,:));
    
    subplot(121);
    hold on;
    h(i) = plot(t,K,'.','MarkerSize',8,'Color',cmap(i,:));
end

subplot(121);
xlim([0 0.5]); ylim([0.3 0.8]);
set(gca,'fontsize',14);
xlabel('$1/\sqrt{t}$ (ms$^{-1/2}$)','Interpreter','latex','FontSize',24);
ylabel('$D(t)/D_0$','Interpreter','latex','FontSize',24);
pbaspect([1 1.5 1]); box on;

subplot(122);
xlim([0 0.5]); ylim([0 2]);
set(gca,'fontsize',14);
xlabel('$1/\sqrt{t}$ (ms$^{-1/2}$)','Interpreter','latex','FontSize',24);
ylabel('$K(t)$','Interpreter','latex','FontSize',24);
pbaspect([1 1.5 1]); box on;
legend(h,{'$\theta=0^\circ$','$\theta=15^\circ$',...
    '$\theta=30^\circ$','$\theta=45^\circ$'},...
    'Interpreter','latex','FontSize',16,'Location','east');

% Plot Fig. 2g-h
mu = [0 0 1];                   % Mean direction: z-axis
% Concentration parameter: no dispersion to high dispersion
kappa = [0.05:0.05:1 1.1:0.1:10 11:30 32:2:100 Inf]; kappa = kappa(end:-1:1);
Ntheta = numel(kappa);          % # polar angle
Nphi = 1;                       % # azimuthal angle

cos2 = zeros(Ntheta,1);         % Mean cosine square of axon directions
Dinf = zeros(Ntheta,2);         % Bulk diffusivity at long time
c = zeros(Ntheta,2);            % Strength of restrictions

simi = sim(1,:);
cmap = colormap('lines');
cmap = cmap([1 5:7],:);
for i = 1:Ntheta
    % Sampling directions based on the Watson distribution
    rng(0);
    if k == 1
        n = repmat([0 0 1],numel(simi)*Nphi,1);
    else
        n = randWatson(numel(simi)*Nphi, mu, kappa(i));
    end
    cos2(k) = mean(n(:,3).^2);
    
    % Calculate displacement cumulants
    dx2i = []; dx4i = [];
    for j = 1:numel(simi)
        list = (j-1)*Nphi+1 : j*Nphi;
        ni = n(list,:);
        [Kj,Dj] = simi(j).data.akc_mom(ni);
        t = simi(j).data.TD;
        dx2j = Dj*2.*t;
        dx2j = mean(dx2j,2);
        dx2i = cat(2,dx2i,dx2j);
    end
    dx2 = sum(fi.*dx2i,2);
    
    % Calculate diffusivity
    D = dx2/2./t;
    
    % Fit the structural disorder model to AD = Dinf + c/sqrt(t)
    tlist = 200:800;
    X = [ones(numel(tlist),1) 1./sqrt(t(tlist))]\D(tlist);
    Dinf(i,1) = X(1);
    c(i,1) = X(2);
    
    % Fit the modified model to AD = Dinf + c/sqrt(t) + c'/t
    X = [ones(numel(tlist),1) 1./sqrt(t(tlist)) 1./t(tlist)]\D(tlist);
    Dinf(i,2) = X(1);
    c(i,2) = X(2);
end

figure('unit','inch','position',[0 0 12 6]);
subplot(121);
hold on;
plot(cos2,Dinf(:,1)/Dinf(1,1),'.','MarkerSize',10);
plot(cos2,Dinf(:,2)/Dinf(1,2),'.','MarkerSize',10);
xlim([0 1]); ylim([0 1]);
box on;
refline(1,0);
xline(1/3);
set(gca,'fontsize',14);
xlabel('$\langle\cos^2\theta\rangle$','Interpreter','latex','FontSize',24);
ylabel('$D_{\infty}/D_{\infty,\,\theta=0}$','Interpreter','latex','FontSize',24);

subplot(122);
hold on;
h1 = plot(cos2,c(:,1)/c(1,1),'.','MarkerSize',10);
h2 = plot(cos2,c(:,2)/c(1,2),'.','MarkerSize',10);
xlim([0 1]); ylim([0 1]);
box on;
refline(1,0);
xline(1/3);
set(gca,'Fontsize',14);
xlabel('$\langle\cos^2\theta\rangle$','Interpreter','latex','FontSize',24);
ylabel('$c/c_{\,\theta=0}$','Interpreter','latex','FontSize',24);
legend([h1 h2],{'uncorrected','corrected'},'Interpreter','latex','FontSize',20);











