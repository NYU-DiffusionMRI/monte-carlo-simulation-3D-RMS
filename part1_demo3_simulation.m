% Part 1: Diffusivity time-dependence along realistic white matter axons
% Demo 3: Simulations (Lee, et al., Commun Biol 2020)
%
% The purpose of this demontration includes:
% 1. Performing simulations of diffusion in four geometries in Part 1 Demo
%    3 (Fig. 2b-f).
%
% Reference: 
% Lee, et al., Communications Biology 2020 (doi:10.1038/s42003-020-1050-x)
% Lee, et al., Brain Structure and Function 2019 
%    (doi:10.1007/s00429-019-01844-6)
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

% Setup directory to this code
root = pwd;

% Load RMS setup tool
addpath(fullfile(root,'analysis'));

%% Create input files of RMS

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

% Load b-table
btable = load(fullfile(root,'hpc_code','btable','btable_5shell_90ndir_b1000.mat'));
bval = btable.bval;
bvec = btable.bvec;

% Setup simulation parameters and microgeometry
% The CUDA code only supports at most 3 compartments by default. To
% increase the number of compartments > 3, please change the Nc_max in line
% 37 of ./hpc_code/lib/main.cu

% The microgeometry has 3 compartments, labeled as 0 (null compartment), 
% 1 (cytoplasm), and 2 (mitochondria) in this example.
% No random walker will be initialized and diffuse in the null compartment.
% The present code supports exchange with the premeation probability from
% compartment a to b: P_{a to b} = min(1,sqrt(D_b/D_a));
sr = setupRMS();
for i = 1:4
    clear X
    X.time_step = 2e-4;                 % Time of each step (ms)
    X.step_num = 5e5;                   % # step
    X.particle_num = 1e3;               % # particle
    X.comp_num = 3;                     % # compartments
    X.voxel_size = file(1).fiber(1).voxel_size;     % Voxel size (micron)
    if i ==1
        % The null compartment must have a diffusivity = 0.
        % The T2 relaxation must be a positive number for all compartments
        X.diffusivity = [0 2 2/15];     % Diffusivity in each compartment (micron2/ms)
        X.T2_relaxation = [80 80 20];   % T2 relaxation time in each compartment (ms)
    else
        X.diffusivity = [0 2 2];        % Diffusivity in each compartment (micron2/ms)
        X.T2_relaxation = [80 80 80];   % T2 relaxation time in each compartment (ms)
    end
    X.bval = bval;                      % b-value (ms/micron2) [Nbvec x 1]
    X.bvec = bvec;                      % gradient directions  [Nbvec x 3]
    for j = 1:numel(file(i).fiber)
        % Setup directory to save the input txt files
        targetj = fullfile(target{i},sprintf('setup%03u',j));
        mkdir(targetj);
        % The microgeometry has 3 compartments, labeled as 0 (null 
        % compartment), 1 (cytoplasm), and 2 (mitochondria) in this example
        fiberj = file(i).fiber(j).shape;
        sr.inputRMS(targetj,fiberj,X);
    end
end

%% Method 1: Compile CUDA C++ code and run simulations

root_cuda = fullfile(root,'hpc_code','lib');
for i = 1:4
    for j = 1:numel(file(i).fiber)
        targetj = fullfile(target{i},sprintf('setup%03u',j));
        cd(targetj);
        copyfile(fullfile(root_cuda,'main.cu'),targetj);
        system('nvcc -arch=sm_70 main.cu -o main_cuda');
        system('./main_cuda');
        system('rm -f main.cu main_cuda');
    end
end

%% Method 2: Create a shell and run the shell in Terminal

root_cuda = fullfile(root,'hpc_code','lib');
fileID = fopen(fullfile(root,'hpc_code','input','job.sh'),'w');
fprintf(fileID,'#!/bin/bash\n');
for i = 1:4
    for j = 1:numel(file(i).fiber)
        targetj = fullfile(target{i},sprintf('setup%03u',j));
        fprintf(fileID,sprintf('cd %s\n',targetj));
        fprintf(fileID,sprintf('cp -a %s .\n',fullfile(root_cuda,'main.cu')));
        fprintf(fileID,'nvcc -arch=sm_70 main.cu -o main_cuda\n');
        fprintf(fileID,'./main_cuda\n');
        fprintf(fileID,'rm -f main.cu main_cuda\n\n');
    end
end
fclose(fileID);

% In terminal, please get to the directory to the job.sh file and run the 
% shell with "sh job.sh"



