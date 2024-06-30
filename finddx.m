function [dxmax,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,thres);
% finddx.m
%
% This function is part of the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
%
% Paper:
% J. Irving and R. Knight (2006), Numerical modeling of ground-penetrating radar
% in 2-D using MATLAB, Computers & Geosciences, 32(9), 1247–1258.
% doi: https://doi.org/10.1016/j.cageo.2005.11.006.
%
% This function finds the maximum spatial discretization that can be used in the
% 2-D FDTD modeling codes TM_model2d.m and TE_model2d.m, such that numerical 
% dispersion is avoided.  Second-order accurate time and fourth-order-accurate 
% spatial derivatives are assumed (i.e., O(2,4)).  Consequently, 5 field points 
% per minimum wavelength are required.
%
% Note:  The dx value obtained with this program is needed to compute the maximum
% time step (dt) that can be used to avoid numerical instability.  However, the
% time vector and source pulse are required in this code to determine the highest
% frequency component in the source pulse.  For this program, make sure to use a fine
% temporal discretization for the source pulse, such that no frequency components
% present in the pulse are aliased.
%
% Syntax:  [dx,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,thres)
%
% where dxmax = maximum spatial discretization possible (m)
%       wlmin = minimum wavelength in the model (m)
%       fmax = maximum frequency contained in source pulse (Hz)
%       epmax = maximum relative dielectric permittivity in grid
%       mumax = maximum relative magnetic permeability in grid
%       srcpulse = source pulse for FDTD simulation
%       t = associated time vector (s)
%       thres = threshold to determine maximum frequency in source pulse
%               (default = 0.02)
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

if nargin==4; thres=0.02; end

% convert relative permittivity and permeability to true values
mu0 = 1.2566370614e-6;
ep0 = 8.8541878176e-12;
epmax = epmax*ep0;
mumax = mumax*mu0;

% compute amplitude spectrum of source pulse and corresponding frequency vector
n = 2^nextpow2(length(srcpulse));
W = abs(fftshift(fft(srcpulse,n)));
W = W./max(W);
fn = 0.5/(t(2)-t(1));
df = 2.*fn/n;
f = -fn:df:fn-df;
W = W(n/2+1:end);
f = f(n/2+1:end);

% determine the maximum allowable spatial disretization
% (5 grid points per minimum wavelength are needed to avoid dispersion)
fmax = f(max(find(W>=thres)));
wlmin = 1/(fmax*sqrt(epmax*mumax));
dxmax = wlmin/5;
