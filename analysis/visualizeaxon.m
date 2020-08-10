classdef visualizeaxon < handle
    
    properties (GetAccess = public, SetAccess = public)
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
    end
    
    methods (Access = public)
        function this = visualizeaxon()
            
        end
        
        function [h,TR] = plotaxon(this,BW)
            [f,v] = isosurface(BW,0);
            TR = triangulation(f,v);
            h = trisurf(TR);
        end
        
        function [N,CENTERS] = histogram(this,X,EDGES)
            CENTERS = EDGES(1:end-1)/2 + EDGES(2:end)/2;
            N = histcounts(X,EDGES);
            N = N/sum(N)./diff(EDGES);     
        end
        
        function fb = bwcalibervar(this,BW)
            [nx,ny,nz] = size(BW);
            fb = zeros(nx,ny,nz,'logical');
            ri = sqrt(sum(sum(BW,1),2)/pi); ri = ri(:);
            xc = nx/2; yc = ny/2;
            [yi,xi] = meshgrid(1:ny,1:nx);
            r2 = (xi-xc).^2 + (yi-yc).^2;
            for k = 1:nz
                fbk = zeros(nx,ny,'logical');
                fbk(r2 <= ri(k)^2) = 1;
                fb(:,:,k) = fbk;
            end
        end
        
        function fb = bwundulation(this,BW)
            [nx,ny,nz] = size(BW);
            fb = zeros(nx,ny,nz,'logical');
            for k = 1:nz
                fbk = zeros(nx,ny,'logical');
                [I,J] = ind2sub([nx,ny],find(BW(:,:,k)));
                fbk(round(mean(I)),round(mean(J))) = 1;
                fb(:,:,k) = fbk;
            end
            ri = ceil(sqrt(nnz(BW)/nz/pi));

            se = strel('sphere',ri);
            zpad = ri;
            fba = repmat(fb(:,:,1),1,1,zpad); 
            fbb = repmat(fb(:,:,end),1,1,zpad);
            fb = cat(3,fba,fb,fbb);
            fb = convn(fb,se.Neighborhood,'same');
            fb = fb(:,:,zpad+1:end-zpad);
            val = unique(fb);
            vol = zeros(numel(val),1);
            for j = 1:numel(val)
                vol(j) = nnz(fb>val(j));
            end
            [~,I] = min(abs(vol-nnz(BW)));
            fb = logical(fb>val(I));
        end
    end
    
        
end