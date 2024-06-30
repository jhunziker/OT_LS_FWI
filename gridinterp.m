function A2 = gridinterp(A,x,z,x2,z2,method)
% gridinterp.m
% 
% This function is part of the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
%
% This function interpolates the electrical property matrix A (having row and column position
% vectors x and z, respectively) to form the matrix A2 (having row and column position vectors
% x2 and z2, respectively).  The interpolation method by default is nearest neighbour ('nearest').  
% Otherwise the method can be specified as 'linear', 'cubic', or 'spline'.
%
% Syntax:  A2 = gridinterp(A,x,z,x2,z2,method)
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

if nargin==5; method = 'nearest'; end

% transpose for interpolation
A = A';

% transform position vectors into matrices for interp2
[x,z] = meshgrid(x,z);
[x2,z2] = meshgrid(x2,z2);

% perform the interpolation
A2 = interp2(x,z,A,x2,z2,method);

% transpose back for output
A2 = A2';
