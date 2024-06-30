function w = OT_nograd(p,tobs,wobs,wpred,dofig)
% Calculate the OT distance based on the algorithm by Malcolm Sambridge et al. (2022). 
% 
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

fs = 12; % Fontsize
lw = 2; % Linewidth

tpred = tobs;

% nugrid = 80; % Discretization of fingerprint window in vertical direction (amplitude)
% nugrid = 20; % Discretization of fingerprint window in vertical direction (amplitude)
% ntgrid = 512; % Discretization of fingerprint window in horizontal direction (time)
% ntgrid = 32; % Discretization of fingerprint window in horizontal direction (time)
lambdav = 0.04; % Distance scale parameter for density field (parameter s in Sambridge et al. (2022))
% grid = [tpred(1), tpred(end), -2.0, 3.5, nugrid,ntgrid]; % t0,t1,u0,u1,nu,ny
grid = [tpred(1), tpred(end), 2.0*min(wpred), 2.0*max(wpred), p.nugrid,p.ntgrid]; % t0,t1,u0,u1,nu,ny

% Build Fingerprint of observed waveform
[dfield_obs,irays_obs,xrays_obs,lrays_obs,pos_obs] = distance_grid(tobs,wobs,grid);

if dofig==1
    figure;
    imagesc(dfield_obs)
    colorbar
    title('Nearest distance field to waveform')
    set(gca,'YDir','normal','Fontsize',fs)
end

pdf_obs = exp(-abs(dfield_obs)/lambdav);
pdf_obs_amp = sum(pdf_obs(:));
pdf_obs = pdf_obs/pdf_obs_amp;
cdf_obs = cumsum(reshape(pdf_obs.',p.nugrid*p.ntgrid,1));

if dofig==1
    figure;
    imagesc(pdf_obs)
    colorbar
    title('Probability density function')
    set(gca,'YDir','normal','Fontsize',fs)
end

% Build Fingerprint of predicted waveform and calculate its derivative
[dfield_pred,irays_pred,xrays_pred,lrays_pred,pos_pred] = distance_grid(tpred,wpred,grid);

if dofig==1
    figure;
    imagesc(dfield_pred)
    colorbar
    title('Nearest distance field to waveform')
    set(gca,'YDir','normal','Fontsize',fs)
end

pdf_pred = exp(-abs(dfield_pred)/lambdav);
pdf_pred_amp = sum(pdf_pred(:));
pdf_pred = pdf_pred/pdf_pred_amp;
cdf_pred = cumsum(reshape(pdf_pred.',p.nugrid*p.ntgrid,1));

if dofig==1
    figure;
    imagesc(pdf_pred)
    colorbar
    title('Probability density function')
    set(gca,'YDir','normal','Fontsize',fs)
end

% Calculate marginals of observed waveform
f0_obs = (sum(pdf_obs,1)).';
f1_obs = sum(pdf_obs,2);
fx0_obs = (pos_obs(1,:,1)).';
fx1_obs = pos_obs(:,1,2);
f0_cdf_obs = cumsum(f0_obs);
f1_cdf_obs = cumsum(f1_obs);

if dofig==1
    figure;
    subplot(211)
    plot(fx0_obs,f0_obs)
    subplot(212)
    plot(fx1_obs,f1_obs)
    title('Marginal probability density function')
    set(gca,'Fontsize',fs)
end

% Calculate marginals of predicted waveform
f0_pred = (sum(pdf_pred,1)).';
f1_pred = sum(pdf_pred,2);
fx0_pred = (pos_pred(1,:,1)).';
fx1_pred = pos_pred(:,1,2);
f0_cdf_pred = cumsum(f0_pred);
f1_cdf_pred = cumsum(f1_pred);

if dofig==1
    figure;
    subplot(211)
    plot(fx0_pred,f0_pred)
    subplot(212)
    plot(fx1_pred,f1_pred)
    title('Marginal probability density function')
    set(gca,'Fontsize',fs)
end

% Calculate Marginal Wasserstein distances and derivatives
[w,dwdpbar0,dwdpbar1,dwdt0] = MargWasserstein(f0_cdf_pred,f0_cdf_obs,fx0_pred,fx0_obs,f1_cdf_pred,f1_cdf_obs,fx1_pred,fx1_obs,pdf_pred,pdf_pred_amp,p.nugrid,p.ntgrid);
