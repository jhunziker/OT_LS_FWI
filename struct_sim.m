function SSIMval = struct_sim(imag,ref)
% Calculate the structural similarity index measure (SSIM). 
% https://en.wikipedia.org/wiki/Structural_similarity_index_measure
% https://ieeexplore.ieee.org/document/1292216
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

mux = mean(ref(:));
muy = mean(imag(:));
temp = cov(ref,imag);
varx = temp(1,1);
vary = temp(2,2);
varxy = temp(1,2);
c1 = 0.0;% 10^(-6);
c2 = 0.0;% 10^(-3);
SSIMval = (2*mux*muy+c1)*(2*varxy+c2)/((mux^2+muy^2+c1)*(varx+vary+c2));