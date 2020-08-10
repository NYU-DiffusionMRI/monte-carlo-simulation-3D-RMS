% Part 1: Diffusivity time-dependence along realistic white matter axons
% Demo 1: Power spectrum (Lee, et al., Commun Biol 2020)
%
% The purpose of this demontration includes:
% 1. Calculating the power spectrum of realistic axonal shapes along axons
%    (Fig. 1).
% 2. Visualizing axon shapes (Fig. 6).
% 3. Calculating the surface area and volume of mitochondria (Supplementary
%    Fig. 1).
%
% Reference: 
% Lee, et al., Communications Biology 2020 (doi:10.1038/s42003-020-1050-x)
% Lee, et al., Brain Structure and Function 2019 
%    (doi:10.1007/s00429-019-01844-6)
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)


% Load axon shape
root = '.';
file = load(fullfile(root,'cell_segmentation','axon_shape.mat'));

% Load visualization tool
addpath(fullfile(root,'analysis'));


% axon.IAS: Intra-axonal space (including cytoplasm and mitochondria) of
%    227 white mater axons segmented from the eletrcon microsocpy data of a 
%    mouse brain.
% axon.mitochondria: Mitochondria of segmented axons
% The axon's main direction is aligned along the z-axis
axon = file.axon;

% voxel_size: The voxel size of the eletron microscopy data.
voxel_size = file.voxel_size;

%% Calculate and normalize axon radius

% Calculate radius of equivalent circle of axon cross sections
stat = struct([]);
for i = 1:size(axon)
    ias = axon(i).IAS;
    ri = sqrt(sum(sum(ias,1),2)/pi)*sqrt(prod(voxel_size(1:2)));
    stat(i,1).ias_r = ri(:);
end
radius = [stat.ias_r];          % Radius

% Normalize radius using mean cross sectional area
[nz,nax] = size(radius);
ias_vol_mean = pi*sum(radius(:).^2)/nax;
radius_norm = zeros(nz,nax);    % Normalized radius
for i = 1:nax
    ri = radius(:,i);
    vi = pi*sum(ri.^2);
    radius_norm(:,i) = ri*sqrt(1/vi*ias_vol_mean);
end

%% 1D power spectrum
% Concatenate distance between beads
bead_pos = [];                  % Bead position
for i = 1:nax
    ri = radius(:,i);
    ri = smooth(ri,10);
    [pks,locs] = findpeaks(ri,'MinPeakDistance',5,'MinPeakWidth',5);
    temp = zeros(nz,1);
    temp(locs) = 1;
    di = diff(locs)*voxel_size(3);
    bead_pos = cat(2,bead_pos,temp);
end

% Concatenate axons of the same length (18 micron)
bead_distance = [];             % Distance between beads
nite = 100;                     % # iterations
rng(0);
for i = 1:nite
    I = randperm(nax);
    pos = bead_pos(:,I);
    dis = diff(find(pos))*voxel_size(3);
    bead_distance = cat(1,bead_distance,dis(:));
end
a_mean = mean(bead_distance);   % Mean of bead distance
a_std = std(bead_distance);     % Standard deviation of bead distance

% 1D power spectrum
L = nax*nz*voxel_size(3);       % Length of concatenated length
power_spectrum_1d = zeros(nax*nz,1);
nite = 1000;                    % # iterations
rng(0);
for i = 1:nite
    I = randperm(nax);
    rho = bead_pos(:,I); rho = rho(:);
    rhok2 = abs(fft(rho)).^2/L;
    power_spectrum_1d = power_spectrum_1d + rhok2;
end
% 1D power spectrum is normalized by multiplying mean(a)
power_spectrum_1d = power_spectrum_1d/nite*a_mean;
kz = (0:numel(power_spectrum_1d)-1)/L;
kz_a_1d = kz*a_mean;          % kz*mean(a)/(2pi)

%% 3D power spectrum
rng(0);
[nz,nax] = size(radius);
L = nax*nz*voxel_size(3);       % Length of concatenated length

power_spectrum_3d = zeros(nax*nz,1);
nite = 1000;                    % # iterations
for i = 1:nite
    I = randperm(nax);
    ri = radius_norm(:,I); ri = ri(:);
    rho = pi*ri.^2*voxel_size(3);
    vol = sum(rho);
    area = vol/L;
    rhok = abs(fft(rho)).^2/vol/area; rhok = rhok(:);
    power_spectrum_3d = power_spectrum_3d + rhok;
end
% 3D power spectrum is normalized by dividing mean(a)
power_spectrum_3d = power_spectrum_3d/nite/a_mean;

kz = (0:numel(power_spectrum_3d)-1)/L;
kz_a_3d = kz*a_mean;          % kz*mean(a)/(2pi)

%% Plot 1D and 3D power spectrum (Fig. 1c-f)
close all
figure('unit','inch','position',[0 0 12 12]);
[nz,nax] = size(radius);
L = nz*voxel_size(3);

% Plot radius variation along 7 selected axons
rng(1);
I = randperm(nax,7);
subplot(221)
z = (1:nz)*voxel_size(3); z = z(:);
h = plot(z,radius(:,I)); set(h,'linewidth',1.5)
xlim([0 L]); ylim([0 1.5])
set(gca,'fontsize',16); axis square
xlabel('z-axis ($\mu$m)','interpreter','latex','fontsize',24);
ylabel('$r$ ($\mu$m)','interpreter','latex','fontsize',24);

% Plot normalized radius variation along 7 selected axons
subplot(222)
h = plot(0.1*(1:180),radius_norm(:,I)); set(h,'linewidth',1.5)
xlim([0 L]); ylim([0 1.5])
set(gca,'fontsize',16); axis square
xlabel('z-axis ($\mu$m)','interpreter','latex','fontsize',24);
ylabel('$\tilde{r}$ ($\mu$m)','interpreter','latex','fontsize',24);

% Plot the histogram of the bead distance
va = visualizeaxon();
edges = linspace(0,25,50);
[N,centers] = va.histogram(bead_distance,edges);
subplot(223)
bar(centers,N,1);
set(gca,'fontsize',16,'xtick',0:5:30,'ytick',0:0.05:0.2)
xlim([0 20]); ylim([0 0.2]); grid on; axis square
xlabel('$a$ ($\mu$m)','interpreter','latex','fontsize',24);
ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',24);
fprintf('Distance between local swelling is: mean %.2f um and std %.2f um. \n',a_mean,a_std);

% Plot 1D power spectrum
subplot(224)
yyaxis right
plot(kz_a_1d,power_spectrum_1d,'linewidth',1); 
set(gca,'xscale','log','yscale','log'); grid on
xlim([1e-2 1e1]); ylim([1e-4 1.001e1]);
set(gca,'fontsize',16,'xtick',[0.01 0.1 1 10]);
xlabel('$k_z\bar{a}/2\pi$','interpreter','latex','fontsize',24);
ylabel('$\Gamma_{\rm pos}(k_z)\cdot \bar{a}$','interpreter','latex','fontsize',24);

% Plot 3D power spectrum. The part of kz < 1/(individual axon length) is
% not shown because the length scale longer than individual axon length is
% unrelated with the bead arrangement.
yyaxis left
La = nz*voxel_size(3);        % Individual axon length
[~,I1] = min(abs(kz_a_3d-a_mean/La));
kz_a_3d(1:I1-1) = nan;
plot(kz_a_3d,power_spectrum_3d,'linewidth',1);
% plateau: The level of power spectrum plateau at low k
plateau = mean(power_spectrum_3d(I1:I1+10));
hr = refline(0,plateau); set(hr,'linestyle','--','linewidth',0.75)
ylim([1e-4 1.001e1]);
set(gca,'xscale','log','yscale','log'); grid on
ylabel('$\Gamma_{\rm 1d}(k_z)/\bar{a}$','interpreter','latex','fontsize',24);
pbaspect([1 1 1]);


%% Visualize axon shapes (Fig. 6a)

figure('unit','inch','position',[0 0 15 5])
j = 0;
rng(8)
va = visualizeaxon();
for i = randsample(numel(axon),6).'
    ias = axon(i).IAS; 
    mit = axon(i).mitochondria;
    
    j = j+1;
    subplot(1,6,j)
    
    h = va.plotaxon(ias);
    set(h,'edgealpha',0,'facealpha',0.6,'facecolor',[0.3 0.3 0.3])
    axis equal off
    hold on
    material dull; lightangle(-45,30)
    
    if nnz(mit)>0
        h = va.plotaxon(mit);
        set(h,'edgealpha',0,'facealpha',0.6,'facecolor',[1 0 0])
        material dull; 
    end
    lightangle(-45,30)
end

%% Relation between mitochondria and axon caliber variation (Fig. 6b and c)

stat = struct([]);
for i = 1:numel(axon)
    ias = axon(i).IAS; 
    mit = axon(i).mitochondria;
    
    ri = sqrt(sum(sum(ias,1),2)/pi)*sqrt(prod(voxel_size(1:2))); 
    ri = ri(:).';
    % Mean radius of IAS for each axon
    stat(i,1).ias_r = mean(ri);
    
    Im = find(sum(sum(mit,1),2));
    Ic = setdiff(1:size(ias,3),Im);
    % Radius of IAS cross sections without the presence of mitochondria
    stat(i,1).ias_cyt_r = ri(Ic);
    
    % Radius of IAS cross sections with the presence of mitochondria
    stat(i,1).ias_mit_r = ri(Im);
    
    % Mitochondria volume for each axon
    stat(i,1).mit_vol = nnz(mit)*prod(voxel_size);
    
    % Axon length for each axon
    stat(i,1).ias_len = size(ias,3)*voxel_size(3);
end

ias_cyt_r = [stat.ias_cyt_r]; ias_cyt_r = ias_cyt_r(:);
ias_mit_r = [stat.ias_mit_r]; ias_mit_r = ias_mit_r(:);
ias_d = 2*[stat.ias_r]; ias_d = ias_d(:);
mit_vol = [stat.mit_vol]; mit_vol = mit_vol(:);
ias_len = [stat.ias_len]; ias_len = ias_len(:);

va = visualizeaxon();
edges = linspace(0,4,51);
[Na,centers] = va.histogram(2*ias_cyt_r,edges);
[Nm,centers] = va.histogram(2*ias_mit_r,edges);

figure('unit','inch','position',[0 0 10 5]);
% Axonal diameter vs mitochondrial volume per unit axonal length
subplot(121)
plot(ias_d,mit_vol./ias_len,'.','markersize',8);
set(gca,'xscale','log','yscale','log');
box on; grid on

I0 = find(mit_vol==0);
dtmp = ias_d; vtmp = mit_vol./ias_len;
dtmp(I0) = []; vtmp(I0) = [];
X = [ones(numel(dtmp),1) log(dtmp)]\log(vtmp);

hold on
fit_d = linspace(0.5,2.2,100);
fit_v = exp(X(1))*fit_d.^X(2);
hr = plot(fit_d,fit_v,'-r');
set(hr,'linewidth',1);
set(gca,'xtick',0.5:0.5:3,'fontsize',12)
xlabel('Inner axonal diameter ($\mu$m)','interpreter','latex','fontsize',16);
ylabel('V(mitocondria)/axonal length ($\mu$m$^2$)','interpreter','latex','fontsize',16);
xlim([0.5 2.2])
pbaspect([1 1 1])

% Axonal diameter of cross sections with and without the presence of
% mitochondria
subplot(122)
hold on;
ha = bar(centers,Na,1); set(ha,'facealpha',0.5);
hm = bar(centers,Nm,1); set(hm,'facealpha',0.5);
legend([ha,hm],{'w/o mitochondria','w/ mitochondria'},'interpreter','latex','fontsize',16)
box on; grid on;
set(gca,'ytick',0:0.5:5,'xtick',0:4,'fontsize',12)
xlim([0 4]); ylim([0 1.5]); pbaspect([1 1 1])
xlabel('Inner axonal diameter ($\mu$m)','interpreter','latex','fontsize',16);
ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
p = ranksum(2*ias_cyt_r,2*ias_mit_r,'method','approximate','tail','right');

%% Histogram of mitochondiral surface area and volume (Supplementary Fig. 1)
% In Supplementary Fig. 1, the volume and surface area were estimated
% before aligning axons to the z-axis. However, these metrics here are
% estimated based on the well-aligned geometry. This will lead to minor
% differences in the histogram.

stat = struct([]);
for i = 1:numel(axon)
    ias = axon(i).IAS; 
    mit = axon(i).mitochondria;
    
    xrange = find(sum(sum(mit,2),3)); xrange = xrange(:);
    yrange = find(sum(sum(mit,1),3)); yrange = yrange(:);
    zrange = find(sum(sum(mit,1),2)); zrange = zrange(:);
    mit = mit(xrange,yrange,zrange);
    
    CC = bwconncomp(mit);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    vol = numPixels*prod(voxel_size);
    % Exclude segmentations of volume smaller than 0.01 micron^3
    vol = vol(vol>0.01);
    
    % Mitochondrial number for each axon
    stat(i,1).mit_num = numel(vol);
    
    % Mitochondrial volume for each mitochondria
    stat(i,1).mit_vol_ind = vol;
    
    % Mitochondrial volume for each axon
    stat(i,1).mit_vol = sum(vol);
    
    % IAS volume for each axon
    stat(i,1).ias_vol = nnz(ias)*prod(voxel_size);
    
    s = regionprops3(mit,'Volume','SurfaceArea');
    vol = s.Volume*prod(voxel_size);
    area = s.SurfaceArea*prod(voxel_size(1:2));
    area = area(vol>0.01); vol = vol(vol>0.01);
    % Mitochondrial surface area for each mitochondria
    stat(i,1).mit_surface_area = area.';
    
    % Ratio of mitochondrial surface area to cytoplasm volume for each axon
    stat(i,1).mit_surface_cyt_volume_ratio = ...
        sum(area)/(nnz(ias)*prod(voxel_size)-sum(vol));
end

mit_num = [stat.mit_num]; mit_num = mit_num(:);
mit_vol = [stat.mit_vol]; mit_vol = mit_vol(:);
ias_vol = [stat.ias_vol]; ias_vol = ias_vol(:);
mit_vol_ind = [stat.mit_vol_ind]; mit_vol_ind = mit_vol_ind(:);

mit_area = [stat.mit_surface_area]; mit_area = mit_area(:);
mit_sv_ratio = [stat.mit_surface_cyt_volume_ratio];
mit_sv_ratio = mit_sv_ratio(:);

va = visualizeaxon();

figure('unit','inch','position',[0,0,15,10]);
% Number of mitochondria per unit IAS volume for all axons
subplot(231);
edges = linspace(0,1,51);
[N,centers] = va.histogram(mit_num./ias_vol,edges);
bar(centers,N,1);
set(gca,'xtick',0:0.25:2,'ytick',0:10,'fontsize',12)
xlabel('\# mitochondria/V(IAS) ($\mu$m$^{-3}$)','interpreter','latex','fontsize',16);
ylabel('PDF ($\mu$m$^3$)','interpreter','latex','fontsize',16)
pbaspect([1,1,1]); grid on
xlim([0 1]); ylim([0 5]);

subplot(232)
% Ratio of mitochondrial surface area to cytoplasm volume for all axons
edges = linspace(0,2,51);
[N,centers] = va.histogram(mit_sv_ratio,edges);
bar(centers,N,1);
set(gca,'xtick',0:0.5:2,'ytick',0:0.5:2,'fontsize',12)
xlabel('S(mitochondria)/V(cytoplasm) ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
ylabel('PDF ($\mu$m)','interpreter','latex','fontsize',16)
pbaspect([1,1,1]); grid on
xlim([0 2]); ylim([0 2]);

subplot(233);
% Ratio of mitochondrial volume to IAS volume for all axons
edges = linspace(0,0.2,51);
[N,centers] = va.histogram(mit_vol./ias_vol,edges);
N = N/sum(N)/mean(diff(edges));
bar(centers,N,1);
set(gca,'xtick',0:0.05:0.2,'ytick',0:5:30,'fontsize',12)
xlabel('V(mitochondria)/V(IAS)','interpreter','latex','fontsize',16);
ylabel('PDF','interpreter','latex','fontsize',16)
pbaspect([1,1,1]); grid on
xlim([0 0.2]); ylim([0 25]);

subplot(234)
% Surface area of all segmented mitochondria
edges = linspace(0,20,51);
[N,centers] = va.histogram(mit_area,edges);
bar(centers,N,1);
set(gca,'xtick',0:5:20,'ytick',0:0.1:0.5,'fontsize',12)
xlabel('$s$(mitochondria) ($\mu$m$^2$)','interpreter','latex','fontsize',16);
ylabel('PDF ($\mu$m$^{-2}$)','interpreter','latex','fontsize',16)
pbaspect([1,1,1]); grid on
xlim([0 20]); ylim([0 0.4]);

subplot(235);
% Mitochondrial volume of all segmented mitochondria
edges = linspace(0,2,51);
[N,centers] = va.histogram(mit_vol_ind,edges);
N = N/sum(N)/mean(diff(edges));
bar(centers,N,1);
set(gca,'xtick',0:0.5:2,'ytick',0:6,'fontsize',12)
xlabel('$v$(mitochondria) ($\mu$m$^3$)','interpreter','latex','fontsize',16);
ylabel('PDF ($\mu$m$^{-3}$)','interpreter','latex','fontsize',16)
pbaspect([1,1,1]); grid on
xlim([0 2]); ylim([0 6]);




