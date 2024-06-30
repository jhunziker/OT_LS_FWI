% Invert crosshole GPR data using the OT-LS-FWI algorithm.
% Paper: J. Hunziker, G. Meles, and N. Linde (2024), Crosshole ground-penetrating radar
% full-waveform inversion by combining optimal-transport and least-squares distances,
% Journal of Applied Geophysics, in submitted.
% Source: https://github.com/jhunziker/OT-LS-FWI 
% 
% This package uses the optimal-transport algorithm of Malcolm Sambridge 
% (Australian National University)
% Paper: 
% M. Sambridge, A. Jackson, and A. P. Valentine (2022), Geophysical inversion and optimal
% transport, Geophysical Journal International, 231(1), 172–198.
% doi: https://doi.org/10.1093/gji/ggac151
% Source code: 
% https://github.com/msambridge/waveform-ot 
%
% This package uses also the 2D GPR crosshole finite difference time domain (FDTD)
% solver of James Irving (University of Lausanne).
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
clear all; close all; clc; tic; 

% Source and receiver positions
srcz = linspace(1,5.5,10).';
srcx = ones(size(srcz))*0.5;
recz = linspace(1,5,41).';
recx = ones(size(recz))*3.5;
nsrc = length(srcz);
nrec = length(recz);

% Model parameters
nx = 41; % Elements in model domain in x-direction
nz = 61; % Elements in model domain in z-direction
x = linspace(0,4,nx); % Coordinate vector in x-direction
z = linspace(0,6,nz); % Coordinate vector in z-direction
epmin = 4.0; % Minimum relative electric permittivity below which no models can be proposed [-]
% Calculate dt based on this value. If a smaller value were proposed, then the FDTD code fails, because dt is too large. 
epmax = 22.0; % Maximum relative electric permittivity above which no models can be proposed [-]
% Calculate dx based on this value. If a larger value were proposed, then the FDTD code fails, because dx is too large. 
% Conductivity
for isrc=1:nsrc
    for irec=1:nrec
        params(isrc,irec).sig = 0.001; % Electric conductivity (assumed known) [S/m]
    end
end

% Starting model in relative electric permittivity [-]
ep = ones(nx,nz)*6;

% Discretization parameters
for isrc=1:nsrc
    for irec=1:nrec
        params(isrc,irec).dxFDTD = 0.04;
        params(isrc,irec).npml = 10; % number of PML boundary layers
        params(isrc,irec).nugrid = 40; % Vertical discretization of fingerprint window (amplitude)
        params(isrc,irec).ntgrid = 128; % Horizontal discretization of fingerprint window (time)
    end
end

% Source and receiver location
for isrc=1:nsrc
    for irec=1:nrec
        params(isrc,irec).srcx = srcx(isrc); % x-coordinate of source position [m]
        params(isrc,irec).srcz = srcz(isrc); % z-coordinate of source position [m]
        params(isrc,irec).recx = recx(irec); % x-coordinate of receiver position [m]
        params(isrc,irec).recz = recz(irec); % z-coordinate of receiver position [m]
        params(isrc,irec).src = 2; % Source type: 1 --> Ex, 2 --> Ez
    end
end

% Parameters for the numerical calculation of the gradient
pmax = 64; % Amount of random points for gradient calculation
pgrow = 1.015; % Factor with which pmax is increased
maxfrac = 0.25; % Maximum fraction of pixels that can be master points
critdist = 0.2; % Distance around the antennas where no point is allowed to be. 
perturb = 0.1; % Perturbation of relative electric permittivity model

% Inversion parameters
maxit = 256; % Total amount of iterations
critfrac = 0.7; % If more than that fraction of traces are less than half a period away from the observed data, the inversion switches from OT to LS. 
srcfrac = 0.5; % Fraction of sources used in mini-batch
recfrac = 0.25; % Fraction of receivers used in mini-batch
start_line_search_log = -5; % Logarithm of lower end of step-length values to be tested
end_line_search_log = 10; % Logarithm of upper end of step-length values to be tested
nsteps = 16; % Number of step-length values to be tested. Use a multiple of CoreNum.
nrefine = 3; % Number of refinement steps for step-length search

% Other parameters
fileind = 'testrun'; % File indicator used in filenames
CoreNum = 16; % Number of cores (parallel computing happens only with the brute-force line-search)
doOT = 1; % Do the first iteration an OT step

%------------------------[END OF PARAMETER DEFINITION]------------------------%

% Seed random number generator to draw points for which the gradient is calculated. 
rng('shuffle')

% Load observed data
load dobs.mat tobs dobs 
wobs = dobs;
clear dobs

% Load standard deviation of noise
load dobs.mat noisestd

% Load source pulse
load dobs.mat srcpulse
t = tobs;
for isrc=1:nsrc
    for irec=1:nrec
        params(isrc,irec).srcpulse = srcpulse;
    end
end
clear srcpulse

% Select sources and receivers for mini-batch
temp = randperm(nsrc);
nsrc_mb = round(srcfrac*nsrc);
srcel_vec = sort(temp(1:nsrc_mb));
temp = randperm(nrec);
nrec_mb = round(recfrac*nrec);
recel_vec = sort(temp(1:nrec_mb));
clear temp

% Calculate synthetic data based on starting model (the data are calculated for all sources and receivers for later plotting)
wpred = forward_multisrc_multirec(x,z,ep,params(1,1).sig,srcx,srcz,recx,recz,params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);

% Save data related to starting model
save([fileind,'_starting_model.mat'],'t','wpred','tobs','wobs','x','z','ep','params');

% Calculate horizontal and vertical roughness of starting model
d1 = diff(ep,[],1);
rough_horz = sum(abs(d1(:)));
d2 = diff(ep,[],2);
rough_vert = sum(abs(d2(:)));

% Calculate distance
w_OT = 0.0;
w_LS = 0.0;
for is=1:nsrc_mb
    for ir=1:nrec_mb
        isrc = srcel_vec(is);
        irec = recel_vec(ir);
        w_OT = w_OT + sum(OTonly_nograd(params(isrc,irec),tobs,squeeze(wobs(:,irec,isrc)).',squeeze(wpred(:,irec,isrc)).',0));
        w_LS = w_LS + LSonly_nograd(squeeze(wobs(:,irec,isrc)),squeeze(wpred(:,irec,isrc)));
    end
end

% Both the OT and LS distance are calculated, but only one of them is used.
if doOT==1
    w = w_OT;
else
    w = w_LS;
end

% Open parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(CoreNum);
end

% Calculate gradient explicitly. 
% Select randomly the elements used in the calculation of the gradient.
[xelvec_bf,zelvec_bf] = get_points_for_grad(pmax,critdist,x,z,srcx,srcz,recx,recz);
dwdmvec = zeros(size(xelvec_bf));
parfor ip=1:length(xelvec_bf) % Loop over master points
    % Perturb one model parameter
    ep_bf = ep;
    ep_bf(xelvec_bf(ip),zelvec_bf(ip)) = ep_bf(xelvec_bf(ip),zelvec_bf(ip)) + perturb;
    % Forward simulation for perturbed model
    wpred_bf = forward_multisrc_multirec(x,z,ep_bf,params(1,1).sig,srcx(srcel_vec),srcz(srcel_vec),recx(recel_vec),recz(recel_vec),params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);
    % Distance measure for perturbed model using sources and receivers in mini-batch
    w_bf = 0.0;
    for is=1:nsrc_mb
        for ir=1:nrec_mb
            isrc = srcel_vec(is);
            irec = recel_vec(ir);
            if doOT==1
                w_bf = w_bf + sum(OTonly_nograd(params(isrc,irec),tobs,squeeze(wobs(:,irec,isrc)).',squeeze(wpred_bf(:,ir,is)).',0));
            else
                w_bf = w_bf + LSonly_nograd(squeeze(wobs(:,irec,isrc)),squeeze(wpred_bf(:,ir,is)));
            end
        end
    end
    % Gradient for one master point
    dwdmvec(ip) = (w_bf-w)/(ep_bf(xelvec_bf(ip),zelvec_bf(ip))-ep(xelvec_bf(ip),zelvec_bf(ip)));
end
% Interpolate the gradient on a regular grid
[xgrid,zgrid] = ndgrid(x,z);
dwdm_mat = griddata(x(xelvec_bf),z(zelvec_bf),dwdmvec,xgrid,zgrid,'cubic');
dwdm = reshape(dwdm_mat,[length(x)*length(z),1]);

% Print current distance measure to screen
if doOT
    fprintf(['OT-Iteration: 0, W = ',num2str(w),'\n'])
else
    fprintf(['LS-Iteration: 0, W = ',num2str(w),'\n'])
end
wall = w;

% Set model-update direction to negative gradient.
d = -dwdm;

% Calculate the weighted root mean square error (WRMSE)
temp = wobs-wpred;
wrmse = sqrt(mean((temp(:)/noisestd).^2));

% Initialize vectors that are filled during inversion iterations
it = 1;
w_vec = zeros(1,maxit+1);
w_vec(1) = wall;
w_vec_OT = zeros(1,maxit+1);
w_vec_OT(1) = w_OT;
w_vec_LS = zeros(1,maxit+1);
w_vec_LS(1) = w_LS;
wrmse_vec = zeros(1,maxit+1);
wrmse_vec(1) = wrmse;
rough_horz_vec = zeros(1,maxit+1);
rough_horz_vec(1) = rough_horz;
rough_vert_vec = zeros(1,maxit+1);
rough_vert_vec(1) = rough_vert;
Tfracvec_out = zeros(1,maxit);
m = reshape(ep,[length(x)*length(z),1]);
% Loop over inversion iterations. 
% The inversion stops when the maximum amount of iterations is reached
% or when the WRMSE is smaller than 1.05. 
while (it<maxit && wrmse>1.05)
    % Search step length in update direction using line search 
    for refit=1:nrefine
        if refit==1
            alpha_start_log = start_line_search_log;
            alpha_end_log = end_line_search_log;
            alpha_vec = logspace(alpha_start_log,alpha_end_log,nsteps);
        else
            % Line search in the area of the best step length found in the previous refinement step
            [~,el] = min(temp_vec);
            elstart = el-1;
            if elstart<1
                elstart=1;
            end
            alpha_start_log = log10(alpha_vec(elstart));
            elend = el+1;
            if elend>length(alpha_vec)
                elend=length(alpha_vec);
            end
            alpha_end_log = log10(alpha_vec(elend));
            alpha_vec = logspace(alpha_start_log,alpha_end_log,nsteps);
        end
        % Try out all step-length values of this refinement step
        temp_mat = zeros(nsrc_mb,nrec_mb,length(alpha_vec));
        parfor ia = 1:length(alpha_vec)
            % Propose model
            mprop = m+alpha_vec(ia)*d;
            mprop(mprop<epmin) = epmin;
            mprop(mprop>epmax) = epmax;
            mpropmat = reshape(mprop,[length(x),length(z)]);
            % Forward simulation using proposed model
            wpred = forward_multisrc_multirec(x,z,mpropmat,params(1,1).sig,srcx(srcel_vec),srcz(srcel_vec),recx(recel_vec),recz(recel_vec),params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);
            % Calculate distance measure
            for is=1:nsrc_mb
                for ir=1:nrec_mb
                    isrc = srcel_vec(is);
                    irec = recel_vec(ir);
                    if doOT==1
                        temp_mat(is,ir,ia) = sum(OTonly_nograd(params(isrc,irec),tobs,squeeze(wobs(:,irec,isrc)).',squeeze(wpred(:,ir,is)).',0));
                    else
                        temp_mat(is,ir,ia) = LSonly_nograd(squeeze(wobs(:,irec,isrc)),squeeze(wpred(:,ir,is)));
                    end
                end
            end
        end
        temp_vec = squeeze(sum(sum(temp_mat,1),2)).';
    end

    % Select the best step-length alpha
    [~,el] = min(temp_vec);
    alpha_best = alpha_vec(el);
   
    % Update model 
    wall_old = wall;
    m = m+alpha_best*d;
    m(m<epmin) = epmin;
    m(m>epmax) = epmax;
    mmat = reshape(m,[length(x),length(z)]);

    % Calculate horizontal and vertical roughness of starting model
    d1 = diff(mmat,[],1);
    rough_horz = sum(abs(d1(:)));
    d2 = diff(mmat,[],2);
    rough_vert = sum(abs(d2(:)));

    % Calculate synthetic data based on updated model (calculated on all sources and receivers for calculation of the relative crosscorrelation coefficient and for writing to disk)
    wpred = forward_multisrc_multirec(x,z,mmat,params(1,1).sig,srcx,srcz,recx,recz,params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);

    % Calculate switching criterion
    Tfrac = relative_crosscorrelation_coefficient(wobs,wpred,tobs,0,0);
    Tfracvec = Tfrac(:);
    Tfrac_counter = 0;
    for isrp=1:length(Tfracvec)
        if Tfracvec(isrp)<0.5
            Tfrac_counter = Tfrac_counter + 1;
        end
    end
    if Tfrac_counter/length(Tfracvec)>critfrac
        % Switch to LS
        doOT = 0;
        % From now on, use all the data at each update. 
        srcfrac = 1.0; % Fraction of sources used in mini-batch
        recfrac = 1.0; % Fraction of receivers used in mini-batch
    end
    fprintf([num2str(Tfrac_counter/length(Tfracvec)*100,'%.2f'),'%% of traces are close enough to do LS.\n'])
    % Store the relative crosscorrelation coefficient fraction for writing to disk. 
    Tfracvec_out(it) = Tfrac_counter/length(Tfracvec); 

    % Select sources and receivers for mini-batch
    temp = randperm(nsrc);
    nsrc_mb = round(srcfrac*nsrc);
    srcel_vec = sort(temp(1:nsrc_mb));
    temp = randperm(nrec);
    nrec_mb = round(recfrac*nrec);
    recel_vec = sort(temp(1:nrec_mb));
    clear temp

    % Calculate distance
    w_OT = 0.0;
    w_LS = 0.0;
    for is=1:nsrc_mb
        for ir=1:nrec_mb
            isrc = srcel_vec(is);
            irec = recel_vec(ir);
            w_OT = w_OT + sum(OTonly_nograd(params(isrc,irec),tobs,squeeze(wobs(:,irec,isrc)).',squeeze(wpred(:,irec,isrc)).',0));
            w_LS = w_LS + LSonly_nograd(squeeze(wobs(:,irec,isrc)),squeeze(wpred(:,irec,isrc)));
        end
    end

    if doOT==1
        w = w_OT;
    else
        w = w_LS;
    end

    % Calculate gradient numerically (brute-force approach)
    % Increase pmax
    pmax = round(pmax*pgrow);
    if pmax>maxfrac*nx*nz
        pmax = round(maxfrac*nx*nz);
    end
    % Select randomly the elements used in the calculation of the gradient.
    [xelvec_bf,zelvec_bf] = get_points_for_grad(pmax,critdist,x,z,srcx,srcz,recx,recz);
    dwdmvec = zeros(size(xelvec_bf));
    parfor ip=1:length(xelvec_bf)
        % Pertrub one model parameter
        mmat_bf = mmat;
        mmat_bf(xelvec_bf(ip),zelvec_bf(ip)) = mmat_bf(xelvec_bf(ip),zelvec_bf(ip)) + perturb;
        % Forward simulation for perturbed model
        wpred_bf = forward_multisrc_multirec(x,z,mmat_bf,params(1,1).sig,srcx(srcel_vec),srcz(srcel_vec),recx(recel_vec),recz(recel_vec),params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);
        % Distance measure for perturbed model using sources and receivers in mini-batch
        w_bf = 0.0;
        for is=1:nsrc_mb
            for ir=1:nrec_mb
                isrc = srcel_vec(is);
                irec = recel_vec(ir);
                if doOT==1
                    w_bf = w_bf + sum(OTonly_nograd(params(isrc,irec),tobs,squeeze(wobs(:,irec,isrc)).',squeeze(wpred_bf(:,ir,is)).',0));
                else
                    w_bf = w_bf + LSonly_nograd(squeeze(wobs(:,irec,isrc)),squeeze(wpred_bf(:,ir,is)));
                end
            end
        end
        % Numerical gradient for one master point
        dwdmvec(ip) = (w_bf-w)/(mmat_bf(xelvec_bf(ip),zelvec_bf(ip))-mmat(xelvec_bf(ip),zelvec_bf(ip)));
    end
    % Interpolate the gradient on a regular grid
    [xgrid,zgrid] = ndgrid(x,z);
    dwdm_mat = griddata(x(xelvec_bf),z(zelvec_bf),dwdmvec,xgrid,zgrid,'cubic');
    dwdm = reshape(dwdm_mat,[length(x)*length(z),1]);

    % Store the current value of the optimal-transport distance
    wall = w;

    % Print current distance measure to screen  
    if doOT==1 
        fprintf(['OT-Iteration: ',num2str(it),', W = ',num2str(w),'\n'])
    else
        fprintf(['LS_Iteration: ',num2str(it),', W = ',num2str(w),'\n'])
    end

    % Calculate the weighted root mean square error (WRMSE)
    temp = wobs-wpred;
    wrmse = sqrt(mean((temp(:)/noisestd).^2));

    % Store W for plotting
    w_vec(it+1) = wall;
    w_vec_OT(it+1) = w_OT; 
    w_vec_LS(it+1) = w_LS;
    wrmse_vec(it+1) = wrmse;
    rough_horz_vec(it+1) = rough_horz;
    rough_vert_vec(it+1) = rough_vert;

    % Set update direction
    d = -dwdm;

    if mod(it,16)==0
        % Save intermediate results
        save([fileind,'_results_it',num2str(it),'.mat'],'t','wpred','tobs','wobs','x','z','m','w_vec','w_vec_OT','w_vec_LS','wrmse_vec','rough_horz_vec','rough_vert_vec','Tfracvec_out','params');
    end

    it = it + 1;
end
% If the inversion ended before the maximum amount of iterations was reached, 
% remove the empty elements of vectors for storing values for each iteration. 
w_vec(it+1:end) = [];
w_vec_OT(it+1:end) = [];
w_vec_LS(it+1:end) = [];
wrmse_vec(it+1:end) = [];
rough_horz_vec(it+1:end) = [];
rough_vert_vec(it+1:end) = [];
Tfracvec_out(it:end) = [];

% Print why algorithm stopped to screen.
if wrmse<=noisestd
    fprintf('Algorithm converged (WRMSE<=std(noise)).\n')
else
    fprintf('Algorithm stopped, because maximum amount of iterations is reached.\n')
end

% Forward problem for final model.
mmat = reshape(m,[length(x),length(z)]);
wpred = forward_multisrc_multirec(x,z,mmat,params(1,1).sig,srcx,srcz,recx,recz,params(1,1).src,t,params(1,1).srcpulse,params(1,1).dxFDTD,params(1,1).npml);

% Save final results
save([fileind,'_results_OT.mat'],'t','wpred','tobs','wobs','x','z','m','w_vec','w_vec_OT','w_vec_LS','wrmse_vec','rough_horz_vec','rough_vert_vec','Tfracvec_out','params');

toc
