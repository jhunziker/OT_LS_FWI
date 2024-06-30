% Load the mat-files that contain the results of the test of the steplength
clear all; close all; clc; tic; 

% load test_steplength.mat
% load test_steplength_withfilt.mat

runname = 'testrun01';
iter = 2;
nrefine = 3; % Amount of refinements
load([runname,'_iter',num2str(iter),'_refine1_alpha_temp_vec.mat']); w_vec = temp_vec;

% The following block is only necessary for the 8pieces scripts: 
% Define eight regions within which all pixels will have the same value. 
% selcube = zeros(length(x),length(z),8);
% selcube( 1: 20, 1: 15,1) = 1.0;
% selcube( 1: 20,16: 30,2) = 1.0;
% selcube( 1: 20,31: 45,3) = 1.0;
% selcube( 1: 20,46:end,4) = 1.0;
% selcube(21:end, 1: 15,5) = 1.0;
% selcube(21:end,16: 30,6) = 1.0;
% selcube(21:end,31: 45,7) = 1.0;
% selcube(21:end,46:end,8) = 1.0;

% selcube = zeros(length(x),length(z),24);
% selcube( 1: 10, 1: 10, 2) = 1.0;
% selcube( 1: 10,11: 20, 1) = 1.0;
% selcube( 1: 10,21: 30, 3) = 1.0;
% selcube( 1: 10,31: 40, 4) = 1.0;
% selcube( 1: 10,41: 50, 5) = 1.0;
% selcube( 1: 10,51:end, 6) = 1.0;
% selcube(11: 20, 1: 10, 7) = 1.0;
% selcube(11: 20,11: 20, 8) = 1.0;
% selcube(11: 20,21: 30, 9) = 1.0;
% selcube(11: 20,31: 40,10) = 1.0;
% selcube(11: 20,41: 50,11) = 1.0;
% selcube(11: 20,51:end,12) = 1.0;
% selcube(21: 30, 1: 10,13) = 1.0;
% selcube(21: 30,11: 20,14) = 1.0;
% selcube(21: 30,21: 30,15) = 1.0;
% selcube(21: 30,31: 40,16) = 1.0;
% selcube(21: 30,41: 50,17) = 1.0;
% selcube(21: 30,51:end,18) = 1.0;
% selcube(31:end, 1: 10,19) = 1.0;
% selcube(31:end,11: 20,20) = 1.0;
% selcube(31:end,21: 30,21) = 1.0;
% selcube(31:end,31: 40,22) = 1.0;
% selcube(31:end,41: 50,23) = 1.0;
% selcube(31:end,51:end,24) = 1.0;

fs = 12; % Fontsize
lw = 2; % Linewidth

figure(1); 
for iref=1:nrefine
    load([runname,'_iter',num2str(iter),'_refine',num2str(iref),'_alpha_temp_vec.mat']); w_vec = temp_vec;
    if iref==1
        loglog(alpha_vec,w_vec,'Linewidth',lw)
    else
        hold on
        loglog(alpha_vec,w_vec,'Linewidth',lw)
        hold off
    end
end
xlabel('Alpha [-]')
ylabel('Optimal-Transport distance [-]')
set(gca,'Fontsize',fs)
print('-dpng',[runname,'_iter',num2str(iter),'_alpha.png'])

if exist('temp_smooth_vec','var')
    figure(2);
    for iref=1:nrefine
        load([runname,'_iter',num2str(iter),'_refine',num2str(iref),'_alpha_temp_vec.mat']); w_vec = temp_smooth_vec;
        if iref==1
            loglog(alpha_vec,w_vec,'Linewidth',lw)
        else
            hold on
            loglog(alpha_vec,w_vec,'Linewidth',lw)
            hold off
        end
    end
    xlabel('Alpha [-]')
    ylabel('Smoothness regularization objective function [-]')
    set(gca,'Fontsize',fs)
    print('-dpng',[runname,'_iter',num2str(iter),'_alpha_smoothness.png'])
end

if (exist('temp_smooth_vec','var') && exist('lambda','var'))
    figure(3);
    for iref=1:nrefine
        load([runname,'_iter',num2str(iter),'_refine',num2str(iref),'_alpha_temp_vec.mat']); w_vec = temp_vec; w_smooth_vec = temp_smooth_vec;
        if iref==1
            loglog(alpha_vec,w_vec+lambda*w_smooth_vec,'Linewidth',lw)
        else
            hold on
            loglog(alpha_vec,w_vec+lambda*w_smooth_vec,'Linewidth',lw)
            hold off
        end
    end
    xlabel('Alpha [-]')
    ylabel('Overall objective function [-]')
    set(gca,'Fontsize',fs)
    print('-dpng',[runname,'_iter',num2str(iter),'_alpha_overall.png'])
end

figure; 
imagesc(x,z,reshape(m,[length(x),length(z)]).')
axis image
colorbar
xlabel('Distance (m)')
ylabel('Depth (m)')
title('Model')
set(gca,'Fontsize',fs)
print('-dpng',[runname,'_iter',num2str(iter),'_startmodel.png'])

% The following line is only necessary for the 8pieces scripts: 
% d = reshape(sum(selcube.*repmat(d,[length(x),length(z),1]),3),[length(x)*length(z),1]);

figure; 
imagesc(x,z,reshape(d,[length(x),length(z)]).')
axis image
colorbar
xlabel('Distance (m)')
ylabel('Depth (m)')
title('Gradient')
set(gca,'Fontsize',fs)
print('-dpng',[runname,'_iter',num2str(iter),'_gradient.png'])

% [~,el] = min(w_vec(1,:)+w_vec(2,:));
% 
% figure(2);
% plot(tobs,wobs,'k')
% hold on
% plot(tobs,squeeze(wpredall(el,:)),'r')
% hold off

% for iel=1:size(wpredall,1)
%     figure(2);
%     plot(tobs,wobs,'k')
%     hold on
%     plot(tobs,squeeze(wpredall(iel,:)),'r')
%     hold off
%     pause(0.5)
% end

toc