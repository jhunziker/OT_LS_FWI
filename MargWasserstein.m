function [wpmarg,dwpmargX,dwpmargY,dwgmarg] = MargWasserstein(f0_cdf_pred,f0_cdf_obs,fx0_pred,fx0_obs,f1_cdf_pred,f1_cdf_obs,fx1_pred,fx1_obs,pdf_pred,pdf_pred_amp,nugrid,ntgrid)
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

dwp = zeros(nugrid,ntgrid);
wp = 0.0;
wpmarg  = zeros(2,1); %  W_p for each marginal
dwgmarg = zeros(2,1); % deriv wrt t grid for each marginal

[wsqpd,dw,dwg] = wasser(f0_cdf_pred,f0_cdf_obs,fx0_pred,fx0_obs);
dwp = dwp+repmat(dw.',[nugrid,1]);
dwgmarg(1,1) = dwg; % Why is this not there for the second marginal?
dwpmargX = repmat(dw.',[nugrid,1]);
wpmarg(1,1) = wsqpd;
wp = wp + wsqpd;

[wsqpd,dw,dwg] = wasser(f1_cdf_pred,f1_cdf_obs,fx1_pred,fx1_obs);
dwp = dwp+repmat(dw,[1,ntgrid]);
dwgmarg(2,1) = dwg; % This line was added by jhunziker and maybe shouldn't be there
dwpmargY = repmat(dw,[1,ntgrid]);
wpmarg(2,1) = wsqpd;
wp = wp + wsqpd;
wpav = wp/2;

dwp = dwp - sum(reshape(dwp.',[nugrid*ntgrid,1]).*reshape(pdf_pred.',[nugrid*ntgrid,1]));
dwp = dwp/pdf_pred_amp;
dwpmargX = dwpmargX - sum(reshape(dwpmargX.',[nugrid*ntgrid,1]).*reshape(pdf_pred.',[nugrid*ntgrid,1])); % calculate derivatives w.r.t. unormalised source PDF amplitudes first marginal
dwpmargY = dwpmargY - sum(reshape(dwpmargY.',[nugrid*ntgrid,1]).*reshape(pdf_pred.',[nugrid*ntgrid,1])); % calculate derivatives w.r.t. unormalised source PDF amplitudes second marginal
dwpmargX = dwpmargX/pdf_pred_amp;
dwpmargY = dwpmargY/pdf_pred_amp;
dwpav = dwp/2;
dwgav = dwg/2;
