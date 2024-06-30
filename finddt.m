function dtmax = finddt(epmin,mumin,dx,dz);
% finddt.m
%
% This function is part of the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
%
% This function finds the maximum time step that can be used in the 2-D
% FDTD modeling codes TM_model2d.m and TE_model2d.m, such that they remain
% numerically stable.  Second-order-accurate time and fourth-order-accurate 
% spatial derivatives are assumed (i.e., O(2,4)).
%
% Syntax: dtmax = finddt(epmin,mumin,dx,dz)
%
% where dtmax = maximum time step for FDTD to be stable
%       epmin = minimum relative dielectric permittivity in grid
%       mumin = minimum relative magnetic permeability in grid
%       dx = spatial discretization in x-direction (m)
%       dz = spatial discretization in z-direction (m)
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

% convert relative permittivity and permeability to true values
mu0 = 1.2566370614e-6;
ep0 = 8.8541878176e-12;
epmin = epmin*ep0;
mumin = mumin*mu0;

% determine maximum allowable time step for numerical stability
dtmax = 6/7*sqrt(epmin*mumin/(1/dx^2 + 1/dz^2));
