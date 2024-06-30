% Plot results of run_OT_LS_FWI.m
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
close all; clear all; clc; tic; 

dostartmod = 1; % Show starting model (1) or not(0)
dostartdata = 1; % Show synthetic data related to starting model (1) or not (0)
filename = 'testrun_results_OT.mat'; % File that contains inversion results
filenamestart = 'testrun_starting_model.mat'; % File that contains starting model
file_truemod = 'dobs.mat'; % model wit one big anomaly in the center
fs = 18; % Fontsize
lw = 2; % Linewidth
ms = 12; % Markersize

% Load data corresponding to starting model. 
if dostartdata==1
    load(filenamestart,'wpred');
    wpred_start = wpred;
    clear wpred
end

% Load results
load(filename)

% Plot data
nrecplot = min([size(wobs,2),5]);
recplotvec = round(linspace(1,size(wobs,2),nrecplot));
for ifig=1:size(wobs,3)
    figure(ifig);
    for iplot=1:nrecplot
        subplot(nrecplot,1,iplot)
        plot(tobs*1e9,wobs(:,recplotvec(iplot),ifig),'k','Linewidth',lw)
        hold on
        plot(t*1e9,wpred(:,recplotvec(iplot),ifig),'r--','Linewidth',lw)
        if dostartdata==1
            plot(t*1e9,wpred_start(:,recplotvec(iplot),ifig),'b--','Linewidth',lw)
        end
        hold off
        if iplot==1
            if dostartdata==1
                legend('Observed','Predicted','Starting model','Location','NorthWest')
            else
                legend('Observed','Predicted','Location','NorthWest')
            end
            title(['Data for final model; Source ',num2str(ifig)])
        end
        xlabel('Time (ns)')
        set(gca,'Fontsize',fs)
    end
end

% Plot OT and LS distance as a function of iterations
figure;
semilogy(linspace(0,length(w_vec)-1,length(w_vec)),w_vec_OT,'b','Linewidth',lw)
hold on
semilogy(linspace(0,length(w_vec)-1,length(w_vec)),w_vec_LS,'r','Linewidth',lw)
hold off
xlim([0,length(w_vec)-1])
xlabel('Iteration')
ylabel('W')
legend('OT distance','LS distance')
set(gca,'Fontsize',fs)

% Plot final model
figure; 
imagesc(x,z,reshape(m,[length(x),length(z)]).')
hold on
for isrc=1:size(params,1)
    plot(params(isrc,1).srcx,params(isrc,1).srcz,'*r','Markersize',ms,'Linewidth',lw)
end
for irec=1:size(params,2)
    plot(params(1,irec).recx,params(1,irec).recz,'vr','Markersize',ms,'Linewidth',lw)
end
hold off
axis image
caxis([4,22])
cax = caxis;
colorbar
xlabel('Distance (m)')
ylabel('Depth (m)')
title('Inverted model')
set(gca,'Fontsize',fs)

% Plot starting model
if dostartmod==1
    load(filenamestart,'ep');

    figure;
    imagesc(x,z,ep.')
    hold on
    for isrc=1:size(params,1)
        plot(params(isrc,1).srcx,params(isrc,1).srcz,'*r','Markersize',ms,'Linewidth',lw)
    end
    for irec=1:size(params,2)
        plot(params(1,irec).recx,params(1,irec).recz,'vr','Markersize',ms,'Linewidth',lw)
    end
    hold off
    axis image
    % caxis([3,15])
    caxis(cax);
    colorbar
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    title('Starting model')
    set(gca,'Fontsize',fs)
end

% Plot true model
if exist('file_truemod','var')==1 % Check if the variable exists
    if exist(file_truemod,'file')==2 % Check if the file exists
        load(file_truemod,'ep');

        figure;
        imagesc(x,z,ep.')
        hold on
        for isrc=1:size(params,1)
            plot(params(isrc,1).srcx,params(isrc,1).srcz,'*r','Markersize',ms,'Linewidth',lw)
        end
        for irec=1:size(params,2)
            plot(params(1,irec).recx,params(1,irec).recz,'vr','Markersize',ms,'Linewidth',lw)
        end
        hold off
        axis image
        % caxis([3,15])
        caxis(cax);
        colorbar
        xlabel('Distance (m)')
        ylabel('Depth (m)')
        title('True model')
        set(gca,'Fontsize',fs)
    end
end

% Plot the NRCCC
figure;
plot(linspace(1,length(Tfracvec_out),length(Tfracvec_out)),Tfracvec_out,'k','Linewidth',lw)
xlim([0,length(w_vec)-1])
xlabel('Iteration')
ylabel('Trace fraction ready for LS')
set(gca,'Fontsize',fs)

% Plot WRMSE
figure;
plot(linspace(0,length(wrmse_vec)-1,length(wrmse_vec)),wrmse_vec,'k','Linewidth',lw)
xlim([0,length(wrmse_vec)-1])
xlabel('Iteration')
ylabel('Weighted Root Mean Square Error')
set(gca,'Fontsize',fs)

% Plot roughness
figure;
plot(linspace(0,length(rough_horz_vec)-1,length(rough_horz_vec)),rough_horz_vec,'k','Linewidth',lw)
hold on
plot(linspace(0,length(rough_vert_vec)-1,length(rough_vert_vec)),rough_vert_vec,'k--','Linewidth',lw)
hold off
xlim([0,length(rough_horz_vec)-1])
xlabel('Iteration')
ylabel('Roughness')
legend('Horizontal Roughness','Vertical Roughness','Location','SouthEast')
set(gca,'Fontsize',fs)

% Calculate structural similarity index mesure (SSIM)
mmat = reshape(m,[length(x),length(z)]);
fprintf(['SSIM: ',num2str(struct_sim(mmat,ep)),'\n'])

toc
