function gather = forward_multisrc_multirec(x,z,ep,sig,srcx,srcz,recx,recz,src,t,srcpulse,dxFDTD,npml)
% Calculate the forward problem using the FDTD code of James Irving. 
%
% This function is based on the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
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

% Create a homogeneous model
nx = length(x);
nz = length(z);
mu = ones(nx,nz);
sig = ones(nx,nz)*sig;

% Spatial grid is kept to store snapshots
outputpar.x=x;
outputpar.z=z;
outputpar.returnfields = 0; % return snapshots (1) or not (0)

% calculate minimum and maximum relative permittivity and permeability
% in the model (to be used in finddx.m and finddt.m) 
epmin = min(min(ep));
epmax = max(max(ep));
mumin = min(min(mu));
mumax = max(max(mu));

% use finddx.m to determine maximum possible spatial field discretization
% (in order to avoid numerical dispersion)
[dxmax,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,0.02);
if dxFDTD>dxmax
    error(['Spatial discretization dx for finite-difference simulation has to be smaller than ',num2str(dxmax),' m.'])
end
dx = dxFDTD;
dz = dx;

% find the maximum possible time step using this dx and dz
% (in order to avoid numerical instability)
dtmax = finddt(epmin,mumin,dx,dz);
dt = t(2)-t(1);
if dt>dtmax
    error(['Temporal discretization dt for finite-difference simulation has to be smaller than ',num2str(dtmax),' s.'])
end

% interpolate electrical property grids to proper spatial discretization
% NOTE:  we MUST use dx/2 here because we're dealing with electrical property matrices
x2 = min(x):dx/2:max(x);
z2 = min(z):dx/2:max(z);
% ep2 = gridinterp(ep,x,z,x2,z2,'cubic');
ep2 = gridinterp(ep,x,z,x2,z2,'linear'); % Changed to avoid values outside epmin and epmax
mu2 = gridinterp(mu,x,z,x2,z2,'cubic');
sig2 = gridinterp(sig,x,z,x2,z2,'cubic');

% pad electrical property matrices for PML absorbing boundaries
[ep3,x3,z3] = padgrid(ep2,x2,z2,2*npml+1);
[mu3,x3,z3] = padgrid(mu2,x2,z2,2*npml+1);
[sig3,x3,z3] = padgrid(sig2,x2,z2,2*npml+1);

% clear unnecessary matrices taking up memory
clear x x2 z z2 ep ep2 mu mu2 sig sig2 

% create source and receiver location matrices (includes type)
% (rows are [x location (m), z location (m), type (src: 1=Ex,2=Ez)])
srctype = src*ones(size(srcz));
rectype = src*ones(size(recz));
srcloc = [srcx srcz srctype];
recloc = [recx recz rectype];

% set some output and plotting parameters
doplot = 0; % Plot during simulations (1) or not (0)
outstep = 1;
plotopt = [doplot 2 50 0.001];

% run the simulation
[gather,tout,srcx,srcz,recx,recz,output] = TE_model2d(ep3,mu3,sig3,x3,z3,srcloc,recloc,srcpulse,t,npml,outstep,plotopt,outputpar);
