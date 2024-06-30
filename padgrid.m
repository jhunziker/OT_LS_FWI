function [A2,x2,z2] = padgrid(A,x,z,n)
% padgrid.m
% 
% This function is part of the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
% 
% This function pads the electrical property matrix A (having row and column position vectors
% x and z, respectively) with n elements around each side to create the matrix A2 (having row
% and column position vectors x2 and z2, respectively).  The properties in the padded regions
% are simply the properties of the original matrix extended outwards.
%
% Syntax:  [A2,x2,z2] = padgrid(A,x,z,n)
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

% determine the new position vectors
dx = x(2)-x(1);
dz = z(2)-z(1);
x2 = (x(1)-n*dx):dx:(x(end)+n*dx+1e-10);
z2 = (z(1)-n*dz):dz:(z(end)+n*dz+1e-10);

% pad the grid
A2 = [repmat(A(:,1),1,n), A, repmat(A(:,end),1,n)];
A2 = [repmat(A2(1,:),n,1); A2; repmat(A2(end,:),n,1)];
