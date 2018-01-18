% Emulator workspace

clc; clear all; close all;

%% Get data
addpath('.\NSF DEMS\Phase 1\');

fprintf('Reading data from .xlsx...\n')
raw_dat = xlsread("fe_results.xlsx");
tdat = Tdat(raw_dat,3); % rescaling inputs, standardizing outputs
fprintf('done.\n\n')

%% Get prediction points
% Get a grid of prediction points.
ms=8;
fprintf('Grid mesh set to length out = %g\n\n',ms)

temp_pred  = unique(tdat.input(:,2));
VF_pred    = linspace(0,1,ms);
thick_pred = linspace(0,1,ms);
dum_var    = [0 .5 1];

pred_pts = allcomb(dum_var,temp_pred,VF_pred,thick_pred);

%% Covariance parameter settings
omega = [0.2 0.82];
rho   = [0.39 0.99];
lambda = 0.005;

%% Run the emulator
em=emulator(tdat.input,tdat.output,pred_pts,omega,rho,lambda,0,false);

%% TBDeleted
% ones_vec = ones(size(pred_pts,1),1);
% %% ??? Changing this
% %pred_pts = [ ones_vec ones_vec pred_pts ; 
% %    2*ones_vec ones_vec pred_pts];
% %% to this
% pred_pts = [ 0*ones_vec pred_pts ; .5 * ones_vec pred_pts ; 
%     ones_vec pred_pts ];
% 
% % Rescale inputs, standardize outputs
% temp_sim  = dat(:,1);
% VF_sim    = dat(:,2);
% thick_sim = dat(:,3);
% temp_sim_min   = min(temp_sim);
% temp_sim_range = range(temp_sim);
% VF_sim_min   = min(VF_sim);
% VF_sim_range = range(VF_sim);
% thick_sim_min   = min(thick_sim);
% thick_sim_range = range(thick_sim);
% 
% temp_sim01 = (temp_sim-min(temp_sim))./range(temp_sim);
% VF_sim01 = (VF_sim-min(VF_sim))./range(VF_sim);
% thick_sim01 = (thick_sim-min(thick_sim))./range(thick_sim);
% 
% defl_sim = dat(:,4);
% rot_sim  = dat(:,5);
% cost_sim = dat(:,6);
% defl_sim_mean = mean(defl_sim);
% defl_sim_sd   = std(defl_sim);
% rot_sim_mean = mean(rot_sim);
% rot_sim_sd   = std(rot_sim);
% cost_sim_mean = mean(cost_sim);
% cost_sim_sd   = std(cost_sim);
% 
% defl_sim_std = (defl_sim - defl_sim_mean)/defl_sim_sd;
% rot_sim_std  = (rot_sim - rot_sim_mean)/rot_sim_sd;
% cost_sim_std = (cost_sim - cost_sim_mean)/cost_sim_sd;
% 
% sim_des     = [temp_sim01 VF_sim01 thick_sim01];
% sim_res     = [defl_sim_std rot_sim_std cost_sim_std];
% 
% %[E,G,d] = meshgrid(Epred,Gpred,dpred);
% 
% %ypred = zeros(length(Epred),length(Gpred),length(dpred));
% % This will store the predicted response
% ypred = zeros(size(pred_pts,1),1);
% 
% % First, take a look at the GP without any data
% fprintf('Finding covariance matrix for prediction grid...\n')
% Sigma = gp_cov(omega,pred_pts(:,1:2),pred_pts(:,1:2),...
%     rho,pred_pts(:,3:end),pred_pts(:,3:end));
% fprintf('done.\n\n')
% ypred = mvnrnd(ypred,Sigma);
% % Take a look
% pred_pts0 = pred_pts(pred_pts(:,1)==0,:);
% pred_pts1 = pred_pts(pred_pts(:,1)==1,:);
% mu0=ypred(pred_pts(:,1)==0);
% mu1=ypred(pred_pts(:,1)==1);
% %plot(pred_pts0(:,2),mu0,'.');
% %plot(pred_pts1(:,2),mu1,'.');
% % Now look at just one slice
% idx_s = all(pred_pts1(:,4)==pred_pts1(1,4),2);
% pred_pts_s = pred_pts1(idx_s,:);
% ypred_s = ypred(idx_s);
% %plot(pred_pts_s,ypred_s,'.');
% 
% % Now add observations.
% % We need the following four covariance matrices.
% ones_vec = ones(size(sim_des,1),1);
% x  = [ 0*ones_vec sim_des; .5* ones_vec sim_des ; 
%     ones_vec sim_des ];
% xp = pred_pts;
% fprintf('Getting three more covariance matrices... ')
% Sigma_xx  = gp_cov(omega,x(:,1:2),x(:,1:2),...
%     rho,x(:,3:end),x(:,3:end));
% fprintf('1,... ')
% Sigma_xpx  = gp_cov(omega,xp(:,1:2),x(:,1:2),...
%     rho,xp(:,3:end),x(:,3:end));
% fprintf('2,... ')
% Sigma_xxp  = gp_cov(omega,x(:,1:2),xp(:,1:2),...
%     rho,x(:,3:end),xp(:,3:end));
% fprintf('3.\n\n')
% Sigma_xpxp  = Sigma;
% WN = eye(size(Sigma_xx)) * 10^(-6);
% % Now get the mean and covariance of the gp.
% inv_Sig_xx = inv(Sigma_xx + WN);
% rcond(inv_Sig_xx)
% rcond(Sigma_xx + WN)
% mu = Sigma_xpx * inv_Sig_xx * sim_res(:);
% cov_gp = Sigma_xpxp - Sigma_xpx * inv_Sig_xx * Sigma_xxp;
