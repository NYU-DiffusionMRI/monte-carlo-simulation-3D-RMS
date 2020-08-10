classdef setupRMS < handle
    
    properties (GetAccess = public, SetAccess = public)
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
    end
    
    methods (Access = public)
        function this = setupRMS()
            
        end
        
        function inputRMS(this,target,fiber,X)
            fid = fopen(fullfile(target,'fiber.txt'),'w');
            for i = 1:size(fiber,1)
                for j = 1:size(fiber,2)
                    fprintf(fid,sprintf('%u ',fiber(i,j,:)));
                end
            end
            fclose(fid);

            % b table
            bval = X.bval; bvec = X.bvec;
            fid = fopen(fullfile(target,'bval.txt'),'w');
            fprintf(fid,sprintf('%.8f\n',bval));
            fclose(fid);

            fid = fopen(fullfile(target,'bvec.txt'),'w');
            for i = 1:size(bvec,1)
                fprintf(fid,sprintf('%.8f\n',bvec(i,:)));
            end
            fclose(fid);

            % Simulation parameters
            dt = X.time_step; TN = X.step_num; NPar = X.particle_num;
            Nbvec = size(bvec,1);
            Ncomp = X.comp_num;
            voxel_size = X.voxel_size;
            [NPix1,NPix2,NPix3] = size(fiber);
            fid = fopen(fullfile(target,'simParamInput.txt'),'w');
            fprintf(fid,sprintf('%e\n',dt));
            fprintf(fid,sprintf('%u\n',TN));
            fprintf(fid,sprintf('%u\n',NPar));
            fprintf(fid,sprintf('%u\n',Nbvec));
            fprintf(fid,sprintf('%u\n',Ncomp));
            fprintf(fid,sprintf('%e\n',voxel_size));
            fprintf(fid,sprintf('%u\n',NPix1));
            fprintf(fid,sprintf('%u\n',NPix2));
            fprintf(fid,sprintf('%u\n',NPix3));
            fclose(fid);
            
            D = X.diffusivity;
            fid = fopen(fullfile(target,'diffusivity.txt'),'w');
            fprintf(fid,sprintf('%.8f\n',D));
            fclose(fid);

            T2 = X.T2_relaxation;
            fid = fopen(fullfile(target,'T2.txt'),'w');
            fprintf(fid,sprintf('%.8f\n',T2));
            fclose(fid);
            
            step_size_max = sqrt(6*max(D)*dt);
            if step_size_max >= voxel_size
                sprintf('Warning, the step size is larger than the voxel size.\n');
            end
        end
    end
    
        
end