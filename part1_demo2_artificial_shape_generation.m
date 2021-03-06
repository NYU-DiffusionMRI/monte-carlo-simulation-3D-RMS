% Part 1: Diffusivity time-dependence along realistic white matter axons
% Demo 2: Artificial shape generation (Lee, et al., Commun Biol 2020)
%
% The purpose of this demontration includes:
% 1. Generating artificial microgeometry based on realistic axonal shape
%    (Fig. 2a).
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

%% I. The original axonal shape (cytoplasm = 1, mitochondria = 2)
target = fullfile(root,'hpc_code','input','I_realistic_axon_mitochondria');
fiber = struct([]);
for i = 1:numel(axon)
    fiber(i).shape = uint8(axon(i).IAS + axon(i).mitochondria);
    fiber(i).voxel_size = voxel_size(1);
end
save(fullfile(target,'fiber.mat'),'fiber');

%% II. The original axonal shape (cytoplasm = 1, mitochondria = 1)
target = fullfile(root,'hpc_code','input','II_realistic_axon');
fiber = struct([]);
for i = 1:numel(axon)
    fiber(i).shape = axon(i).IAS;
    fiber(i).voxel_size = voxel_size(1);
end
save(fullfile(target,'fiber.mat'),'fiber');

%% III. Fibers with only caliber variation (no undulation)
target = fullfile(root,'hpc_code','input','III_caliber_variation');
va = visualizeaxon();
fiber = struct([]);
tic;
for i = 1:numel(axon)
    ias = axon(i).IAS;
    fiber(i).shape = va.bwcalibervar(ias);
    fiber(i).voxel_size = voxel_size(1);
end
toc;
save(fullfile(target,'fiber.mat'),'fiber');

%% IV. Fibers with only undulation (no caliber variation)
target = fullfile(root,'hpc_code','input','IV_undulation');
va = visualizeaxon();
fiber = struct([]);
tic;
for i = 1:numel(axon)
    ias = axon(i).IAS;
    fiber(i).shape = va.bwundulation(ias);
    fiber(i).voxel_size = voxel_size(1);
end
toc;
save(fullfile(target,'fiber.mat'),'fiber');

%% Plot fibers (Fig. 2a)
close all
target = cell(4,1);
target{1} = fullfile(root,'hpc_code','input','I_realistic_axon_mitochondria');
target{2} = fullfile(root,'hpc_code','input','II_realistic_axon');
target{3} = fullfile(root,'hpc_code','input','III_caliber_variation');
target{4} = fullfile(root,'hpc_code','input','IV_undulation');

file = struct([]);
for i = 1:4
    filei = load(fullfile(target{i},'fiber.mat'));
    file(i,1).fiber = filei.fiber;
end

figure('unit','inch','position',[0 0 3 10]);
va = visualizeaxon();
xrange = [35 65]; yrange = [35 65]; zrange = [0 180];
j = 0;
rng(0);
for i = randsample(227,3).'
    j = j+1;
    
    % I. The original axonal shape (cytoplasm = 1, mitochondria = 2)
    hs = subplot(3,4,(j-1)*4+1);
    set(hs,'position',[0 (3-j)/3 0.95/4 0.95/3]);
    hold on;
    ha = va.plotaxon(file(1).fiber(i).shape==1);
    hm = va.plotaxon(file(1).fiber(i).shape==2);
    set(ha,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    set(hm,'edgealpha',0,'facecolor','r')
    axis equal off
    xlim(xrange); ylim(yrange); zlim(zrange); pbaspect([range(xrange) range(yrange) range(zrange)]);
    view(3)
    material dull
    light('Position',[1 0 0],'style','local')
    
    % II. The original axonal shape (cytoplasm = 1, mitochondria = 1)
    hs = subplot(3,4,(j-1)*4+2);
    set(hs,'position',[1/4 (3-j)/3 0.95/4 0.95/3]);
    ha = va.plotaxon(file(2).fiber(i).shape);
    set(ha,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    axis equal off
    xlim(xrange); ylim(yrange); zlim(zrange); pbaspect([range(xrange) range(yrange) range(zrange)]);
    view(3)
    material dull
    light('Position',[1 0 0],'style','local')
    
    % III. Fibers with only caliber variation (no undulation)
    hs = subplot(3,4,(j-1)*4+3);
    set(hs,'position',[2/4 (3-j)/3 0.95/4 0.95/3]);
    ha = va.plotaxon(file(3).fiber(i).shape);
    set(ha,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    axis equal off
    xlim(xrange); ylim(yrange); zlim(zrange); pbaspect([range(xrange) range(yrange) range(zrange)]);
    view(3)
    material dull
    light('Position',[1 0 0],'style','local')
    
    % IV. Fibers with only undulation (no caliber variation)
    hs = subplot(3,4,(j-1)*4+4);
    set(hs,'position',[3/4 (3-j)/3 0.95/4 0.95/3]);
    ha = va.plotaxon(file(4).fiber(i).shape);
    set(ha,'edgealpha',0,'facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    axis equal off
    xlim(xrange); ylim(yrange); zlim(zrange); pbaspect([range(xrange) range(yrange) range(zrange)]);
    view(3)
    material dull
    light('Position',[1 0 0],'style','local')
    
end
