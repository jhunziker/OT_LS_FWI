function [gather,tout,srcx,srcz,recx,recz,output] = TE_model2d(ep,mu,sig,xprop,zprop,srcloc,recloc,srcpulse,t,npml,outstep,plotopt,outputpar)
% TE_model2d.m
%
% This function is part of the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
%  
% This is a 2-D, TE-mode, FDTD modeling program for crosshole GPR and vertical radar profiling.  
% See TE_run_example.m for an example of how to use this code.
%
% Features:
% - convolutional PML absorbing boundaries (Roden and Gedney, 2000)
% - calculations in PML boundary regions have been separated from main modeling region for speed
% - second-order accurate time derivatives, fourth-order-accurate spatial derivatives (O(2,4))
%
% Syntax:  [gather,tout,srcx,srcz,recx,recz] = TE_model2d(ep,mu,sig,xprop,zprop,srcloc,recloc,srcpulse,t,npml,outstep,plotopt)
%
% where gather = output data cube (common source gathers are stacked along the third index)
%       tout = output time vector corresponding to first index in gather matrix (s)
%       srcx = vector containing actual source x locations after discretization (m)
%       srcz = vector containing actual source z locations after discretization (m)
%       recx = vector containing actual receiver x locations after discretization (m)
%       recz = vector containing actual receiver z locations after discretization (m)
%       output = structure containing the electromagnetic fields if they
%       need to be returned for the last source position. 
%       ep,mu,sig = electrical property matrices  (rows = x, columns = z;  need odd # of rows and columns; 
%             already padded for PML boundaries; ep and mu are relative to free space; sig units are S/m)
%       xprop,zprop = position vectors corresponding to electrical property matrices (m)
%       srcloc = matrix containing source locations and type (rows are [x (m), z (m), type (1=Ex,2=Ez])
%       recloc = matrix containing receiver locations (format same as above)
%       srcpulse = source pulse vector (length must be equal to the number of iterations)
%       t = time vector corresponding to srcpulse (used for determining time step and number of iterations) (s)
%       npml = number of PML absorbing boundary cells
%       outstep = write a sample to the output matrix every so many iterations (default=1)
%       plotopt = plot Ex or Ez wavefield during simulation?  
%           (vector = [{0=no, 1=yes}, (1=Ex, 2=Ez) {output every # of iterations}, {colorbar threshold}])
%           (default = [1 2 50 0.05])
%       returnfields = Flag to indicate if the electromagnetic fields should be returned or not {0=no, 1=yes}
%
% Notes:
% - all matrices in this code need to be transposed before plotting, as the first index always refers 
%   to the horizontal direction, and the second index to the vertical direction
% - locations in space are referred to using x and z
% - indices in electric and magnetic field matrices are referred to using (i,j)  [i=horizontal,j=vertical]
% - indices in electrical property matrices are referred to using (k,l)  [k=horizontal,l=vertical]
%
% by James Irving
% July 2005
%
% This file is part of the OT-LS-FWI package. OT-LS-FWI is free software: 
% you can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. 
% OT-LS-FWI is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details. You should have received a copy of the GNU General Public 
% License along with OT-LS-FWI. If not, see <https://www.gnu.org/licenses/>.  

% ------------------------------------------------------------------------------------------------
% SET DEFAULTS, AND TEST A FEW THINGS TO MAKE SURE ALL THE INPUTS ARE OK
% ------------------------------------------------------------------------------------------------

if nargin==10; outstep=1; plotopt=[1 2 50 0.05]; end
if nargin==11; plotopt=[1 2 50 0.05]; end

if size(mu)~=size(ep) | size(sig)~=size(ep); disp('ep, mu, and sig matrices must be the same size'); return; end
if [length(xprop),length(zprop)]~=size(ep); disp('xprop and zprop are inconsistent with ep, mu, and sig'); return; end
if mod(size(ep,1),2)~=1 | mod(size(ep,2),2)~=1; disp('ep, mu, and sig must have an odd # of rows and columns'); return; end
if size(srcloc,2)~=3 | size(recloc,2)~=3; disp('srcloc and recloc matrices must have 3 columns'); return; end
if max(srcloc(:,1))>max(xprop) | min(srcloc(:,1))<min(xprop) | max(srcloc(:,2))>max(zprop)...
        | min(srcloc(:,2))<min(zprop); disp('source vector out of range of modeling grid'); return; end 
if max(recloc(:,1))>max(xprop) | min(recloc(:,1))<min(xprop) | max(recloc(:,2))>max(zprop)...
        | min(recloc(:,2))<min(zprop); disp('receiver vector out of range of modeling grid'); return; end
if length(srcpulse)~=length(t); disp('srcpulse and t vectors must have same # of points'); return; end
if npml>=length(xprop)/2 | npml>=length(zprop)/2; disp('too many PML boundary layers for grid'); return; end
if length(plotopt)~=4; disp('plotopt must be a 4 component vector'); return; end


% ------------------------------------------------------------------------------------------------
% DETERMINE SOME INITIAL PARAMETERS
% ------------------------------------------------------------------------------------------------

% determine true permittivity and permeability matrices from supplied relative ones
ep0 = 8.8541878176e-12;             % dielectric permittivity of free space 
mu0 = 1.2566370614e-6;              % magnetic permeability of free space
ep = ep*ep0;                        % true permittivity matrix            
mu = mu*mu0;                        % true permeability matrix

% determine number of field nodes and discretization interval
nx = (length(xprop)+1)/2;           % maximum number of field nodes in the x-direction
nz = (length(zprop)+1)/2;           % maximum number of field nodes in the z-direction
dx = 2*(xprop(2)-xprop(1));         % electric and magnetic field spatial discretization in x (m)
dz = 2*(zprop(2)-zprop(1));         % electric and magnetic field spatial discretization in z (m)

% x and z position vectors corresponding to Ex, Ez, and Hy field matrices
% (these field matrices are staggered in both space and time, and thus have different coordinates)
xEx = xprop(2):dx:xprop(end-1);
zEx = zprop(1):dz:zprop(end);
xEz = xprop(1):dx:xprop(end);
zEz = zprop(2):dz:zprop(end-1);
xHy = xEx;
zHy = zEz;

% determine source and receiver (i,j) indices in Ex or Ez field matrices,
% and true coordinates of sources and receivers in numerical model (after discretization)
nsrc = size(srcloc,1);                              % number of sources
nrec = size(recloc,1);                              % number of receivers
for s=1:nsrc;
    if srcloc(s,3)==1                               % in the case of an Ex source...
        [temp,srci(s)] = min(abs(xEx-srcloc(s,1))); % source x index in Ex field matrix
        [temp,srcj(s)] = min(abs(zEx-srcloc(s,2))); % source z index in Ex field matrix
        srcx(s) = xEx(srci(s));                     % true source x location
        srcz(s) = zEx(srcj(s));                     % true source z location
    elseif srcloc(s,3)==2                           % in the case of an Ez source...
        [temp,srci(s)] = min(abs(xEz-srcloc(s,1))); % source x index in Ez field matrix
        [temp,srcj(s)] = min(abs(zEz-srcloc(s,2))); % source z index in Ez field matrix
        srcx(s) = xEz(srci(s));                     % true source x location
        srcz(s) = zEz(srcj(s));                     % true source z location
    end    
end
for r=1:nrec;                                   
    if recloc(r,3)==1                               % in the case of an Ex receiver
        [temp,reci(r)] = min(abs(xEx-recloc(r,1))); % receiver x index in Ex field matrix
        [temp,recj(r)] = min(abs(zEx-recloc(r,2))); % receiver z index in Ex field matrix
        recx(r) = xEx(reci(r));                     % true receiver x location
        recz(r) = zEx(recj(r));                     % true receiver z location
    elseif recloc(r,3)==2                           % in the case of an Ez receiver
        [temp,reci(r)] = min(abs(xEz-recloc(r,1))); % receiver x index in Ez field matrix
        [temp,recj(r)] = min(abs(zEz-recloc(r,2))); % receiver z index in Ez field matrix
        recx(r) = xEz(reci(r));                     % true receiver x location
        recz(r) = zEz(recj(r));                     % true receiver z location
    end    
end

% determine time stepping parameters from supplied time vector
dt = t(2)-t(1);                                 % temporal discretization
numit = length(t);                              % number of iterations


% ------------------------------------------------------------------------------------------------
% COMPUTE FDTD UPDATE COEFFICIENTS FOR ENTIRE SIMULATION GRID
% note:  these matrices are twice as large as the field component matrices
% (i.e., the same size as the electrical property matrices)
% ------------------------------------------------------------------------------------------------

% disp('Determining update coefficients for simulation region...')

% set the basic PML parameters, keeping the following in mind...   
% - maximum sigma_x and sigma_z vary in heterogeneous media to damp waves most effiently
% - Kmax = 1 is the original PML of Berenger (1994); Kmax > 1 can be used to damp evanescent waves
% - keep alpha = 0 except for highly elongated domains (see Roden and Gedney (2000))
m = 4;                                          % PML Exponent (should be between 3 and 4)
Kxmax = 5;                                      % maximum value for PML K_x parameter (must be >=1)
Kzmax = 5;                                      % maximum value for PML K_z parameter (must be >=1)
sigxmax = (m+1)./(150*pi*sqrt(ep./ep0)*dx);     % maximum value for PML sigma_x parameter
sigzmax = (m+1)./(150*pi*sqrt(ep./ep0)*dz);     % maximum value for PML sigma_z parameter
alpha = 0;                                      % alpha parameter for PML (CFS)

% indices corresponding to edges of PML regions in electrical property grids
kpmlLout = 1;                           % x index for outside of PML region on left-hand side
kpmlLin = 2*npml+2;                     % x index for inside of PML region on left-hand side
kpmlRin = length(xprop)-(2*npml+2)+1;   % x index for inside of PML region on right-hand side
kpmlRout = length(xprop);               % x index for outside of PML region on right-hand side
lpmlTout = 1;                           % z index for outside of PML region at the top
lpmlTin = 2*npml+2;                     % z index for inside of PML region at the top
lpmlBin = length(zprop)-(2*npml+2)+1;   % z index for inside of PML region at the bottom
lpmlBout = length(zprop);               % z index for outside of PML region at the bottom

% determine the ratio between the distance into the PML and the PML thickness in x and z directions
% done for each point in electrical property grids;  non-PML regions are set to zero
xdel = zeros(length(xprop),length(zprop));                      % initialize x direction matrix
k = kpmlLout:kpmlLin;  k = k(:);                                % left-hand PML layer
xdel(k,:) = repmat(((kpmlLin-k)./(2*npml)),1,length(zprop));
k = kpmlRin:kpmlRout; k = k(:);                                 % right-hand PML layer
xdel(k,:) = repmat(((k-kpmlRin)./(2*npml)),1,length(zprop));
zdel = zeros(length(xprop),length(zprop));                      % initialize z direction matrix
l = lpmlTout:lpmlTin;                                           % top PML layer
zdel(:,l) = repmat(((lpmlTin-l)./(2*npml)),length(xprop),1);
l = lpmlBin:lpmlBout;                                           % bottom PML layer
zdel(:,l) = repmat(((l-lpmlBin)./(2*npml)),length(xprop),1);

% determine PML parameters at each point in the simulation grid 
% (scaled to increase from the inside to the outside of the PML region)
% (interior non-PML nodes have sigx=sigz=0, Kx=Kz=1)
sigx = sigxmax.*xdel.^m;
sigz = sigzmax.*zdel.^m;
Kx = 1 + (Kxmax-1)*xdel.^m;
Kz = 1 + (Kzmax-1)*zdel.^m;

% determine FDTD update coefficients
Ca  = (1-dt*sig./(2*ep))./(1+dt*sig./(2*ep));
Cbx  = (dt./ep)./((1+dt*sig./(2*ep))*24*dx.*Kx);
Cbz  = (dt./ep)./((1+dt*sig./(2*ep))*24*dz.*Kz);
Cc = (dt./ep)./(1+dt*sig./(2*ep));
Dbx = (dt./(mu.*Kx*24*dx));
Dbz = (dt./(mu.*Kz*24*dz));
Dc = dt./mu;
Bx = exp(-(sigx./Kx + alpha)*(dt/ep0));
Bz = exp(-(sigz./Kz + alpha)*(dt/ep0));
Ax = (sigx./(sigx.*Kx + Kx.^2*alpha + 1e-20).*(Bx-1))./(24*dx);
Az = (sigz./(sigz.*Kz + Kz.^2*alpha + 1e-20).*(Bz-1))./(24*dz);

% clear unnecessary PML variables as they take up lots of memory
clear sigxmax sigzmax xdel zdel Kx Kz sigx sigz


% ------------------------------------------------------------------------------------------------
% RUN THE FDTD SIMULATION
% ------------------------------------------------------------------------------------------------

% disp('Beginning FDTD simulation...')

% initialize gather matrix where data will be stored
gather = zeros(fix((numit-1)/outstep)+1,nrec,nsrc);

% loop over number of sources
for s=1:nsrc

    % zero all field matrices
    Hy = zeros(nx-1,nz-1);          % Hy component of magnetic field
    Ex = zeros(nx-1,nz);            % Ex component of electric field
    Ez = zeros(nx,nz-1);            % Ez component of electric field
    Hydiffx = zeros(nx,nz-1);       % difference for dHy/dx
    Hydiffz = zeros(nx-1,nz);       % difference for dHy/dz
    Exdiffz = zeros(nx-1,nz-1);     % difference for dEx/dz
    Ezdiffx = zeros(nx-1,nz-1);     % difference for dEz/dx
    PHyx = zeros(nx-1,nz-1);        % psi_Hyx (for PML)
    PHyz = zeros(nx-1,nz-1);        % psi_Hyz (for PML)
    PEx = zeros(nx-1,nz);           % psi_Ex (for PML)
    PEz = zeros(nx,nz-1);           % psi_Ez (for PML)
    if outputpar.returnfields==1
        [xHygrid,zHygrid] = meshgrid(xHy,zHy);
        [xExgrid,zExgrid] = meshgrid(xEx,zEx);
        [xEzgrid,zEzgrid] = meshgrid(xEz,zEz);
        [output.xgrid,output.zgrid] = meshgrid(outputpar.x,outputpar.z);
        output.Hy = zeros(length(outputpar.x),length(outputpar.z),numit);
        output.Ex = zeros(length(outputpar.x),length(outputpar.z),numit);
        output.Ez = zeros(length(outputpar.x),length(outputpar.z),numit);
        output.t = t;
    else
        output = 0;
    end
    
    % time stepping loop
    for it=1:numit
        
        % update Hy component...
        
        % determine indices for entire, PML, and interior regions in Hy and property grids
        i = 2:nx-2;  j = 2:nz-2;                % indices for all components in Hy matrix to update
        k = 2*i;  l = 2*j;                      % corresponding indices in property grids
        kp = k((k<=kpmlLin | k>=kpmlRin));      % corresponding property indices in PML region
        lp = l((l<=lpmlTin | l>=lpmlBin));
        ki = k((k>kpmlLin & k<kpmlRin));        % corresponding property indices in interior (non-PML) region
        li = l((l>lpmlTin & l<lpmlBin));
        ip = kp./2;  jp = lp./2;                % Hy indices in PML region
        ii = ki./2;  ji = li./2;                % Hy indices in interior (non-PML) region

        % update to be applied to the whole Hy grid
        Exdiffz(i,j) = -Ex(i,j+2) + 27*Ex(i,j+1) - 27*Ex(i,j) + Ex(i,j-1);
        Ezdiffx(i,j) = -Ez(i+2,j) + 27*Ez(i+1,j) - 27*Ez(i,j) + Ez(i-1,j);
        Hy(i,j) = Hy(i,j) + Dbz(k,l).*Exdiffz(i,j) - Dbx(k,l).*Ezdiffx(i,j);
        
        % update to be applied only to the PML region
        PHyx(ip,j) = Bx(kp,l).*PHyx(ip,j) + Ax(kp,l).*Ezdiffx(ip,j);
        PHyx(ii,jp) = Bx(ki,lp).*PHyx(ii,jp) + Ax(ki,lp).*Ezdiffx(ii,jp);
        PHyz(ip,j) = Bz(kp,l).*PHyz(ip,j) + Az(kp,l).*Exdiffz(ip,j);
        PHyz(ii,jp) = Bz(ki,lp).*PHyz(ii,jp) + Az(ki,lp).*Exdiffz(ii,jp);
        Hy(ip,j) = Hy(ip,j) + Dc(kp,l).*(PHyz(ip,j) - PHyx(ip,j));
        Hy(ii,jp) = Hy(ii,jp) + Dc(ki,lp).*(PHyz(ii,jp) - PHyx(ii,jp));
        
        
        % update Ex component...
        
        % determine indices for entire, PML, and interior regions in Ex and property grids
        i = 2:nx-2;  j = 3:nz-2;                % indices for all components in Ex matrix to update
        k = 2*i;  l = 2*j-1;                    % corresponding indices in property grids
        kp = k((k<=kpmlLin | k>=kpmlRin));      % corresponding property indices in PML region
        lp = l((l<=lpmlTin | l>=lpmlBin));
        ki = k((k>kpmlLin & k<kpmlRin));        % corresponding property indices in interior (non-PML) region
        li = l((l>lpmlTin & l<lpmlBin));    
        ip = kp./2;  jp = (lp+1)./2;            % Ex indices in PML region
        ii = ki./2;  ji = (li+1)./2;            % Ex indices in interior (non-PML) region
        
        % update to be applied to the whole Ex grid
        Hydiffz(i,j) = -Hy(i,j+1) + 27*Hy(i,j) - 27*Hy(i,j-1) + Hy(i,j-2);
	    Ex(i,j) = Ca(k,l).*Ex(i,j) + Cbz(k,l).*Hydiffz(i,j);
        
        % update to be applied only to the PML region
        PEx(ip,j) = Bz(kp,l).*PEx(ip,j) + Az(kp,l).*Hydiffz(ip,j);
        PEx(ii,jp) = Bz(ki,lp).*PEx(ii,jp) + Az(ki,lp).*Hydiffz(ii,jp);
        Ex(ip,j) = Ex(ip,j) + Cc(kp,l).*PEx(ip,j);
        Ex(ii,jp) = Ex(ii,jp) + Cc(ki,lp).*PEx(ii,jp);
        
    
        % update Ez component...

        % determine indices for entire, PML, and interior regions in Ez and property grids
        i = 3:nx-2;  j = 2:nz-2;                % indices for all components in Ez matrix to update
        k = 2*i-1;  l = 2*j;                    % corresponding indices in property grids
        kp = k((k<=kpmlLin | k>=kpmlRin));      % corresponding property indices in PML region
        lp = l((l<=lpmlTin | l>=lpmlBin));
        ki = k((k>kpmlLin & k<kpmlRin));        % corresponding property indices in interior (non-PML) region
        li = l((l>lpmlTin & l<lpmlBin));
        ip = (kp+1)./2;  jp = lp./2;            % Ez indices in PML region
        ii = (ki+1)./2;  ji = li./2;            % Ez indices in interior (non-PML) region

        % update to be applied to the whole Ez grid
        Hydiffx(i,j) = -Hy(i+1,j) + 27*Hy(i,j) - 27*Hy(i-1,j) + Hy(i-2,j);
	    Ez(i,j) = Ca(k,l).*Ez(i,j) - Cbx(k,l).*Hydiffx(i,j);
        
        % update to be applied only to the PML region
        PEz(ip,j) = Bx(kp,l).*PEz(ip,j) + Ax(kp,l).*Hydiffx(ip,j);
        PEz(ii,jp) = Bx(ki,lp).*PEz(ii,jp) + Ax(ki,lp).*Hydiffx(ii,jp);
        Ez(ip,j) = Ez(ip,j) - Cc(kp,l).*PEz(ip,j);
        Ez(ii,jp) = Ez(ii,jp) + Cc(ki,lp).*PEz(ii,jp);
        
        

        % add source pulse to Ex or Ez at source location
        % (emulates infinitesimal Ex or Ez directed line source with current = srcpulse)
        i = srci(s); j = srcj(s);
        if srcloc(s,3)==1
            Ex(i,j) = Ex(i,j) + srcpulse(it);
        elseif srcloc(s,3)==2
            Ez(i,j) = Ez(i,j) + srcpulse(it);
        end
    
        
        % plot the Ex or Ez wavefield if necessary
        if plotopt(1)==1;   
            if mod(it-1,plotopt(3))==0
                disp(['Source ',num2str(s),'/',num2str(nsrc),', Iteration ',num2str(it),'/',num2str(numit),...
                        ':  t = ',num2str(t(it)*1e9),' ns'])
                if plotopt(2)==1
                    figure(1); imagesc(xEx,zEx,Ex'); axis image
                    title(['Source ',num2str(s),'/',num2str(nsrc),', Iteration ',num2str(it),'/',num2str(numit),...
                            ':  Ex wavefield at t = ',num2str(t(it)*1e9),' ns']);
                    xlabel('Position (m)');  ylabel('Depth (m)');
                    caxis([-plotopt(4) plotopt(4)]);
                    pause(0.01);
                elseif plotopt(2)==2
                    figure(1); imagesc(xEz,zEz,Ez'); axis image
                    title(['Source ',num2str(s),'/',num2str(nsrc),', Iteration ',num2str(it),'/',num2str(numit),...
                            ':  Ez wavefield at t = ',num2str(t(it)*1e9),' ns']);
                    xlabel('Position (m)');  ylabel('Depth (m)');
                    caxis([-plotopt(4) plotopt(4)]);
                    pause(0.01);
                end
            end
        end
    
        % record the results in gather matrix if necessary
        if mod(it-1,outstep)==0
            tout((it-1)/outstep+1) = t(it);
            for r=1:nrec
                if recloc(r,3)==1
                    gather((it-1)/outstep+1,r,s) = Ex(reci(r),recj(r));
                elseif recloc(r,3)==2
                    gather((it-1)/outstep+1,r,s) = Ez(reci(r),recj(r));
                end
            end
        end
        
        % save wavefield
        if outputpar.returnfields==1
            output.Hy(:,:,it) = (interp2(xHygrid,zHygrid,Hy(1:end-1,1:end-1)',output.xgrid,output.zgrid))';
            output.Ex(:,:,it) = (interp2(xExgrid,zExgrid,Ex(1:end-1,1:end-1)',output.xgrid,output.zgrid))';
            output.Ez(:,:,it) = (interp2(xEzgrid,zEzgrid,Ez(1:end-1,1:end-1)',output.xgrid,output.zgrid))';
        end

    end
end
