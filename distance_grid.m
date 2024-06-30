function [dfield,irays,xrays,lrays,pos] = distance_grid(tsig,wsig,grid)
% This function is based on the optimal-transport algorithm of Malcolm Sambridge 
% (Australian National University)
% Paper: 
% M. Sambridge, A. Jackson, and A. P. Valentine (2022), Geophysical inversion and optimal
% transport, Geophysical Journal International, 231(1), 172–198.
% doi: https://doi.org/10.1093/gji/ggac151
% Source code: 
% https://github.com/msambridge/waveform-ot 
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

nugrid = grid(5);
ntgrid = grid(6); 

tlimn = [(tsig(1)-grid(1))/(grid(2)-grid(1)),(tsig(end)-grid(1))/(grid(2)-grid(1))]; % Waveform time range in normalized coordinates
ulimn = [0.0,1.0];
delgrid = [(ulimn(2)-ulimn(1))/nugrid,(tlimn(2)-tlimn(1))/ntgrid];
p = [tsig;wsig].'; % Coordinates of waveform points in unnormalised coordinates
pn = [(tsig-grid(1))/(grid(2)-grid(1));(wsig-grid(3))/(grid(4)-grid(3))].'; % Coordinates of waveform points in normalised coordinates
x0(1,:,:) = pn(1:end-1,:); % Coordinates of waveform segment start points in normalised coordinates
delta_n = pn(2:end,:)-pn(1:end-1,:); % Vectors along waveform segments in normalised coordinates
lsq_n = sum(delta_n.^2,2); % Waveform segment squared lengths in normalised coordinates

nump = nugrid*ntgrid; 
[Xn,Yn] = meshgrid(linspace(tlimn(1),tlimn(2),ntgrid),linspace(ulimn(1),ulimn(2),nugrid));
Xnt = Xn';
Ynt = Yn';
temp(:,1,:) = [Xnt(:),Ynt(:)];
b = temp-x0;
clear Xnt Ynt temp
delta_n_temp(1,:,:) = delta_n;
lam(:,:,1) = sum(b.*delta_n_temp,3);
lsq_n_temp(1,:,:) = lsq_n;
lam = lam./lsq_n_temp;
clear lsq_n_temp
lam(lam>1.0) = 1.0;
lam(lam<0.0) = 0.0;
lam_temp(:,:,1) = lam;
ds = b-delta_n_temp.*lam_temp;
clear delta_n_temp lam_temp
dsq = sum(ds.^2,3);
[~,irays] = min(dsq.');
lrays = (lam(sub2ind(size(lam),linspace(1,nugrid*ntgrid,nugrid*ntgrid),irays))).';
xrays = squeeze(x0(1,irays,:))+[lrays,lrays].*delta_n(irays,:);
d = (sqrt(dsq(sub2ind(size(dsq),linspace(1,nugrid*ntgrid,nugrid*ntgrid),irays)))).';
dfield = reshape(d,[ntgrid,nugrid]).';
pos(:,:,1) = Xn;
pos(:,:,2) = Yn;
