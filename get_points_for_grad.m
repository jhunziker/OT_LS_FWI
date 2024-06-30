function [xelvec_bf,zelvec_bf] = get_points_for_grad(pmax,critdist,x,z,srcx,srcz,recx,recz)
% Create two vectors containing the element numbers of random points inside the model domain. 
% The gradient is calculated numerically on these points. 
% 
% Input: 
% - pmax: Amount of random points
% - critdist: Distance around antennas in meters within which no points are allowed to be
% - x, z: Coordinates in x- and z-direction in meters
% - srcx, srcz, recx, recz: Source and receiver coordinates in meters
% 
% Output: 
% - xelvec_bf, zelvec_bf: Vectors containing the element number in x- and z-direction of
%   the master points used for the calculation of the gradient. 
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

nx = length(x);
nz = length(z);
% The gradient needs to be calculated at the four corners of the model. 
xelvec_bf = [1,1,nx,nx];
zelvec_bf = [1,nz,1,nz];
% The remaining points are chosen randomly
pcount = length(xelvec_bf);
while pcount<pmax
    redraw = 0;
    xtemp = randi([1,nx],1);
    ztemp = randi([1,nz],1);
    % Check if the point has already been drawn
    for ip=1:length(xelvec_bf)
        if (xtemp==xelvec_bf(ip) && ztemp==zelvec_bf(ip))
            redraw=1;
            break
        end
    end
    if redraw==1
        continue
    end
    % Check if point is close to sources
    for ip=1:length(srcx)
        tempdist = sqrt((srcx(ip)-x(xtemp))^2 + (srcz(ip)-z(ztemp))^2);
        if tempdist<critdist
            redraw=1;
            break
        end
    end
    if redraw==1
        continue
    end
    % Check if point is close to antennas
    for ip=1:length(recx)
        tempdist = sqrt((recx(ip)-x(xtemp))^2 + (recz(ip)-z(ztemp))^2);
        if tempdist<critdist
            redraw=1;
            break
        end
    end
    if redraw==1
        continue
    end
    pcount = pcount+1;
    xelvec_bf(pcount) = xtemp;
    zelvec_bf(pcount) = ztemp;
end
