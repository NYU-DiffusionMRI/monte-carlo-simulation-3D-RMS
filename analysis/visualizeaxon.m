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
    end
    
        
end