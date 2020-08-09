classdef undulation < handle
    properties (Constant = true, Access = protected)
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
    end
    
    properties (GetAccess = private, SetAccess = protected)
        w2_ind = [1 1; 1 2; 2 2];
        w2_cnt = [1 2 1];
        w4_ind = [1 1 1 1;1 1 1 2; 1 1 2 2; 1 2 2 2; 2 2 2 2];
        w4_cnt = [1 2 3 2 1];
    end
    
    methods (Access = public)
        function this = undulation(varargin)
            
        end
        
        function cm = randundul(this,seed,nz,Lz,dx,dy,nstep,span)
            rng(seed)
            xt = dx*cumsum(2*rand(1,nstep)-1);
            yt = dy*cumsum(2*rand(1,nstep)-1);
            xt = imgaussfilt(xt,span); yt = imgaussfilt(yt,span);
            xt = xt-mean(xt); yt = yt-mean(yt);
            xt = imresize(xt(:),[nz,1],'bilinear');
            yt = imresize(yt(:),[nz,1],'bilinear');
            zt = linspace(0,Lz,nz); zt = zt(:);
            
            cm = [xt(:),yt(:),zt(:)];
            u = normr(cm(end,:)-cm(1,:));
            cmr = this.rotateaxis(cm,u,mean(cm,1),'forward');
            zq = zt(:);
            xq = interp1(cmr(:,3),cmr(:,1),zq,'pchip');
            yq = interp1(cmr(:,3),cmr(:,2),zq,'pchip');
            xq = xq-mean(xq); yq = yq-mean(yq);
            cm = [xq,yq,zq];
        end
        
        function cm = lissajous(this,nz,Lz,w,la,phi)
            z = linspace(0,Lz,nz);
            x = w(1)*sin(2*pi/la(1)*z+phi);
            y = w(2)*sin(2*pi/la(2)*z);
            cm = [x(:) y(:) z(:)];
        end
        
        function pb = randbeadpos(this,seed,nz,Lz,abar,astd)
            rng(seed)
            res = Lz/nz;
            pos = zeros(nz*2,1);
            xb = 1;
            while xb <= 2*nz
                randnum = 0;
                while randnum < 1
                    randnum = abar/res + astd/res*randn;
                end
                xb = xb + round(randnum);
                if xb <= 2*nz
                    pos(round(xb)) = 1;
                elseif xb > 2*nz
                    break
                end
            end
            randnum = round(max(1,nz*rand+1));
            pb = pos(randnum:randnum+nz-1);
        end
        
        function rb = randbeadrad(this,nz,Lz,pb,lbar,rcsa,cv)
            if cv==0
                rb = rcsa*ones(nz,1);
            else
                res = Lz/nz;
                sigma = lbar/2.355;
                N = ceil(sigma/res)*6;
                alpha = (N-1)/(sigma/res)/2;
                bead_profile = gausswin(N,alpha);
                rng(0)
                temp = conv(repmat(pb,3,1),bead_profile,'same');
                temp = temp(end/3+1:end/3*2);
                cv_dr = std(temp)/mean(temp);
                c2 = (cv_dr/cv-1);
                sc = sqrt(cv^2*(c2+1)^2 + 1 + 2*c2 + c2^2);
                dr = rcsa/sc*temp/mean(temp(:));
                r0 = c2*mean(dr);
                rb = dr + r0;
            end
        end
        
        function BW = undulrad2bw(this,nz,Lz,cm,rb)
            res = Lz/nz;
            cm = cm/res; rb = rb/res;
            xrange = max(rb) + max(abs(cm(:,1))); xrange = ceil(2*xrange)+4;
            yrange = max(rb) + max(abs(cm(:,2))); yrange = ceil(2*yrange)+4;
            zrange = nz;

            xi = (cm(:,1) + xrange/2);
            yi = (cm(:,2) + yrange/2);
%             zi = 1:nz; zi = zi(:);
            
            BW = false(xrange,yrange,zrange);
            for i = 1:numel(rb)
                % convolution with 2d disc
                for ii = ceil(xi(i)-rb(i)):(ceil(xi(i)+rb(i))+1)
                    if (ii>xrange)||(ii<1)
                        continue;
                    end
                    for jj = ceil(yi(i)-rb(i)):(ceil(yi(i)+rb(i))+1)
                        if (jj>yrange)||(jj<1)
                            continue;
                        end
                        if this.inside_circle(ii,jj,xi(i),yi(i),rb(i))
                            BW(ii,jj,i) = true;
                        end
                    end
                end
                
%                 % convolution with 3d sphere
%                 for ii = ceil(xi(i)-rb(i)):ceil(xi(i)+rb(i))
%                     if (ii>xrange)||(ii<1)
%                         continue;
%                     end
%                     for jj = ceil(yi(i)-rb(i)):ceil(yi(i)+rb(i))
%                         if (jj>yrange)||(jj<1)
%                             continue;
%                         end
%                         for kk = ceil(zi(i)-rb(i)):ceil(zi(i)+rb(i))
%                             if kk>zrange
%                                 tk = kk-zrange;
%                             elseif kk<1
%                                 tk = kk+zrange;
%                             else
%                                 tk = kk;
%                             end
%                             
%                             if this.inside_sphere(ii,jj,kk,xi(i),yi(i),zi(i),rb(i))
%                                 BW(ii,jj,tk) = true;
%                             end
%                         end
%                     end
%                 end
            end
        end
        
        function BW = undulse2bw(this,nz,Lz,cm,se)
            res = Lz/nz;
            cm = cm/res; rb = ceil(length(se)/2);
            xrange = max(rb) + max(abs(cm(:,1))); xrange = ceil(2*xrange)+4;
            yrange = max(rb) + max(abs(cm(:,2))); yrange = ceil(2*yrange)+4;
            zrange = nz;

            xi = (cm(:,1) + xrange/2);
            yi = (cm(:,2) + yrange/2);
            
            BW = false(xrange,yrange,zrange);
            for i = 1:nz
                % convolution with 2d structural element
                BWi = zeros(xrange,yrange,'single');
                BWi(ceil(xi(i)),ceil(yi(i))) = 1;
                BW(:,:,i) = conv2(BWi,se,'same')>0;
            end
        end
        
        function cm = skeleton(this,BW,vox,span)
            cm = this.centermass(BW,vox);
            cm(:,1) = imgaussfilt(cm(:,1),span);
            cm(:,2) = imgaussfilt(cm(:,2),span);
        end
            
        function [w2,w4] = gausssamp(this,cm,TD,D0)
            nt = numel(TD);
            w2 = zeros(nt,3);
            w4 = zeros(nt,5);
            
            ds = sqrt(sum(diff(cm,1,1).^2,2));
            wx = cm(1:end-1,1)/2+cm(2:end,1)/2;
            wy = cm(1:end-1,2)/2+cm(2:end,2)/2;
            wx = wx-mean(wx);
            wy = wy-mean(wy);
            wxx = (wx-wx.').^2;
            wxy = (wx-wx.').*(wy-wy.');
            wyy = (wy-wy.').^2;

            wxxxx = (wx-wx.').^4;
            wxxxy = (wx-wx.').^3 .* (wy-wy.');
            wxxyy = (wx-wx.').^2 .* (wy-wy.').^2;
            wxyyy = (wx-wx.') .* (wy-wy.').^3;
            wyyyy = (wy-wy.').^4;
            
            s = cumsum(ds);
            ss2 = (s-s.').^2;
            dsds = ds*ds.';
            for i = 1:nt
                pg = exp(-ss2/4/D0/TD(i));
                temp = dsds.*pg;
                pgnorm = sum(temp(:));
                temp = dsds.*wxx.*pg/pgnorm;
                w2(i,1) = sum(temp(:));
                temp = dsds.*wxy.*pg/pgnorm;
                w2(i,2) = sum(temp(:));
                temp = dsds.*wyy.*pg/pgnorm;
                w2(i,3) = sum(temp(:));
                
                temp = dsds.*wxxxx.*pg/pgnorm;
                w4(i,1) = sum(temp(:));
                temp = dsds.*wxxxy.*pg/pgnorm;
                w4(i,2) = sum(temp(:));
                temp = dsds.*wxxyy.*pg/pgnorm;
                w4(i,3) = sum(temp(:));
                temp = dsds.*wxyyy.*pg/pgnorm;
                w4(i,4) = sum(temp(:));
                temp = dsds.*wyyyy.*pg/pgnorm;
                w4(i,5) = sum(temp(:));
            end
        end
        
        function z2 = gausssamp_z(this,cm,TD,D0)
            nt = numel(TD);
            z2 = zeros(nt,1);
            
            ds = sqrt(sum(diff(cm,1,1).^2,2));
            z = cm(1:end-1,3)/2+cm(2:end,3)/2;
            zz2 = (z-z.').^2;
            s = cumsum(ds);
            ss2 = (s-s.').^2;
            dsds = ds*ds.';
            for i = 1:nt
                pg = exp(-ss2/4/D0/TD(i));
                temp1 = dsds.*pg;
                temp2 = dsds.*zz2.*pg;
                temp2 = temp2/sum(temp1(:));
%                 temp1 = temp1(floor(end/2-end/20)+1:floor(end/2+end/20),floor(end/2-end/20)+1:floor(end/2+end/20));
%                 temp2 = temp2(floor(end/2-end/20)+1:floor(end/2+end/20),floor(end/2-end/20)+1:floor(end/2+end/20))/sum(temp1(:));
                z2(i) = sum(temp2(:));
            end
        end
        
        function wk = spectrum(this,cm,N)
            wx = cm(1:end-1,1)/2+cm(2:end,1)/2;
            wy = cm(1:end-1,2)/2+cm(2:end,2)/2;
            wx = wx-mean(wx);
            wy = wy-mean(wy);
            wkx = 2*ifft([wx; wx(end:-1:1)]);
            wky = 2*ifft([wy; wy(end:-1:1)]);
            wkx = real(wkx);
            wky = real(wky);
            wk = sqrt(wkx.^2+wky.^2);
            wk = wk(2:N+1);
        end
        
        function [K,D] = akc2d(this,w2,w4,TD)
            phi = (0:179).'*pi/180;
            ni = [cos(phi) sin(phi)];
            ni2 = prod(reshape(ni(:,this.w2_ind),[],3,2),3);
            ni4 = prod(reshape(ni(:,this.w4_ind),[],5,4),3);
            w2i = this.w2_cnt.*w2;
            w4i = this.w4_cnt.*w4;
            w2i = w2i*ni2.';
            w4i = w4i*ni4.';
            wxy2 = mean(w2i,2);
            wxy4 = mean(w4i,2);
            D = wxy2/2./TD(:);
            K = wxy4./wxy2.^2-3;
        end
        
        function Dinst = instdiff(this,D,TD)
            TD = TD(:);
            wxy2 = 4*D.*TD;
            t = TD(1:end-1)/2 + TD(2:end)/2;
            Dinst = 1/4*diff(wxy2,1,1)/mean(diff(t));
            Dinst = interp1(t,Dinst,TD,'spline');
        end
        
        function [Domega,omega] = dispdiff(this,Dinst,TD)
            dt = mean(diff(TD));
%             nt = numel(TD); 
%             x = (1:nt).'; cen = (nt+1)/2;
%             ds = abs(x-cen);
%             dapo = sinc(ds/nt*dt);
            Dw = fft(Dinst)*dt;
%             Dw = Dw./dapo;
            Dw = Dw(ceil(end/2)+1:end);
            omega = (0:numel(Dw)-1).'/max(TD);
            Domega = real(-1j*omega.*Dw);
        end
    end
    
    methods (Static)
        function [xr,M] = rotateaxis(x,u,cm,direction)
            if nargin < 4, direction = 'forward'; end
            v = -[-u(2), u(1), 0];
            s = sqrt(sum(v.^2));
            c = u(3);
            V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            
            R = eye(4);
            R(1:3, 1:3) = eye(3) + V + V*V * (1-c)/s^2;
            if strcmpi(direction,'backward'), R = R'; end
            
            T1 = eye(4); T1(4, 1:3) = cm - [1 1 1]; T1 = T1';
            T2 = eye(4); T2(4, 1:3) = -cm + [1 1 1]; T2 = T2';
            
            M = T1*R*T2;
            if ~isempty(x)
                x = [x.'; ones(1,size(x,1))];
                xr = M*x;
                xr = xr(1:3,:).';
            else
                xr = [];
            end
        end
        
        function inside = inside_sphere(ii,jj,kk,xc,yc,zc,r)
            inside = ( (ii-1/2-xc).^2+(jj-1/2-yc).^2+(kk-1/2-zc).^2 ) <= (r.^2);
        end
        
        function inside = inside_circle(ii,jj,xc,yc,r)
            di = [0 0; 0 1; 1 0; 1 1];
            ti = [1 0; 0 1; -1 0; 0 -1];
            inside = 0;
            for i = 1:size(di,1)
                inside = inside + (( (ii-di(i,1)-xc).^2+(jj-di(i,2)-yc).^2 ) < (r.^2));
                inside = inside + ( (floor(xc+ti(i,1)*r)==(ii-1)) .* (floor(yc+ti(i,2)*r)==(jj-1)) );
            end
            inside = inside>0;
        end
        
        function cm = centermass(BW,vox)
            if numel(vox)==1
                vox = vox*[1 1 1];
            end
            [nx,ny,~] = size(BW);
            zlist = find(sum(sum(BW,1),2));
            cm = zeros(numel(zlist),3);
            for i = 1:numel(zlist)
                fiberi = BW(:,:,zlist(i));
                [I,J] = ind2sub([nx,ny],find(fiberi));
                cm(i,:) = [mean(I)*vox(1),mean(J)*vox(2),zlist(i)*vox(3)];
            end
        end
        
        function savebw(filename,BW)
            save(filename,'BW');
        end
    end
end