classdef readRMS < handle
    
    properties (GetAccess = public, SetAccess = public)
        sig0;
        sig;
        TD;
        dtime; Tstep; NPar; res;
        dx2t; dx4t;
        bval; bvec;
    end
    
    properties (GetAccess = private, SetAccess = private)
        D_cnt = [1 2 2 1 2 1];
        W_cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];
        D_ind = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
        W_ind = [1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3;
                1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3;
                2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
    end
    
    methods (Access = public)
        function this = readRMS(root)
            if ~isempty(root)
                this.readJob(root);
            end
        end
        
        function this = readJob(this,root)
            % Load data
            param = load(fullfile(root,'sim_para.txt'));
            this.dtime = param(1);
            this.Tstep = param(2);
            this.NPar = param(3);
            this.res = param(4);
            this.TD = load(fullfile(root,'diff_time.txt'));
            
            this.bval = load(fullfile(root,'bval.txt'));
            this.bvec = load(fullfile(root,'bvec.txt'));
            this.bvec = reshape(this.bvec,3,[]).';
            
            this.sig0 = load(fullfile(root,'sig0.txt'));
            this.dx2t = load(fullfile(root,'dx2.txt')).*this.D_cnt./this.sig0;
            this.dx4t = load(fullfile(root,'dx4.txt')).*this.W_cnt./this.sig0;
            
            sigRe = load(fullfile(root,'sigRe.txt'));
            this.sig = abs(sigRe)./this.sig0;
        end
        
        function dt = dki_shell(this,sig,bvec,bval)
            Nt = size(sig,1);
            bvecu = unique(bvec,'row');
            n2 = this.ndir2(bvecu);
            n4 = this.ndir4(bvecu);
            
            nbvec = size(bvecu,1);
            nbval = numel(unique(bval));
            dt = zeros(Nt,21);
            for i = 1:Nt
                sigi = sig(i,:);
                Di = zeros(nbvec,1);
                Ki = zeros(nbvec,1);
                for j = 1:nbvec
                    Ij = ismember(bvec,bvecu(j,:),'rows');
                    sigj = sigi(Ij); sigj = sigj(:);
                    bvalj = bval(Ij); bvalj = bvalj(:);
                    A = [-bvalj 1/6*bvalj.^(2:nbval)];
                    X = A\log(sigj + eps);
                    Di(j) = X(1);
                    Ki(j) = X(2)/X(1)^2;
                end
                dx2g = Di*2.*this.TD(i);
                dx4g = (Ki+3).*dx2g.^2;
                dt(i,1:6) = ( n2\dx2g ).';
                dt(i,7:21) = ( n4\dx4g ).';
%                 this.dx2tShl(i,1:6) = ( b(:,1:6)\dx2g ).';
%                 this.dx4tShl(i,7:21) = ( b(:,7:end)\dx4g ).';
            end
        end
        
        function dt = dki_lls(this,sig,bvec,bval)
            % Prepare b table
            grad = [bvec bval];
            normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
            grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
            
            b = -(grad(:,ones(1,6)*4).*grad(:,this.D_ind(1:6,1)).*grad(:,this.D_ind(1:6,2)))*diag(this.D_cnt);
            
            bsqd6 = grad(:, 4).^2 * ones(1,15)/6;
            b = [b, (bsqd6 .* prod(reshape(grad(:,this.W_ind),[],15,4),3))*diag(this.W_cnt)];
            
            % Linear least square
            dtTmp = b\log(sig + eps).';
            dt = DKI.W2K(dtTmp);
%             this.dt_lls = DKI.W2K(dtTmp);
%             
%             % Calculate DTI and DKI metric
%             eigval = DKI.eig(this.dt_lls);
%             this.adSig = DKI.ad(eigval).';
%             this.rdSig = DKI.rd(eigval).';
%             this.mdSig = DKI.md(eigval).';
%             this.akSig = DKI.ak(this.dt_lls).';
%             this.rkSig = DKI.rk(this.dt_lls).';
%             this.mkSig = DKI.mk(this.dt_lls).';
        end
                
        function [K,D] = akc_mom(this,n)
            n2 = this.ndir2(n);
            n4 = this.ndir4(n);
            x2 = this.dx2t*n2.';
            x4 = this.dx4t*n4.';
            D = x2/2./this.TD;
            K = x4./x2.^2-3;
        end
        
        function [K,D] = akc_shell(this,dt,n)
            n2 = this.ndir2(n);
            n4 = this.ndir4(n);
            x2 = dt(:,1:6)*n2.';
            x4 = dt(:,7:21)*n4.';
            D = x2/2./this.TD;
            K = x4./x2.^2-3;
        end
        
        function [K,D] = akc_lls(this,dt,n)
            [K,D] = DKI.akc(dt,n);
            K = K.'; D = D.';
        end
        
        function n4i = ndir4(this,ni)
            n4i = prod(reshape(ni(:,this.W_ind),[],15,4),3);
        end
        
        function n2i = ndir2(this,ni)
            n2i = prod(reshape(ni(:,this.D_ind),[],6,2),3);
        end
    end
    
        
end