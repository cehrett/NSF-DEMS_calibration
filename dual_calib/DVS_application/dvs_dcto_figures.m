% Dynamic Vibration System DCTO figures
% Here I obtain figures for results from dvs_dcto_workspace.m
% 2020-09-25
%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'dual_calib']);
addpath([dpath,'dual_calib\DVS_application']);
addpath([dpath,'dual_calib\DVS_application\data']);
addpath([dpath,'Example']);

% Change dir
cd(dpath);

%% Posterior distributions of inputs
clc ; clearvars -except dpath ; close all ;

% Face alpha for histograms
fa = 0.6;
% Color of true/optimal value line
% lcol = [218 165 32]/255 ; 
lcol='k';

% Load results
loadloc = [dpath,'dual_calib\DVS_application\data\',...
        '\2020-09-30_dvs_dcto_results'];
load(loadloc,'results');

% Define inputs mins and ranges 
theta1 = 6.2e10;
xmin = results.settings.min_x;
xrange = results.settings.range_x;
t1min = results.settings.min_t1;
t1range = results.settings.range_t1;
t2min = results.settings.min_t2;
t2range = results.settings.range_t2;
burn_in = results.settings.burn_in;

% Help function
fillunder = @(x,y,color,falpha) ...
    fill([x(1) x x(end) fliplr(x) 0],...
        [0 y 0 0*y 0],color,'EdgeColor','none','FaceAlpha',falpha);

% First, get prior and posterior theta1
f = figure('pos',[10 10 650 250]);
subplot(1,2,1);
% Plot prior
falpha=.1;
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta1
% histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
%     'EdgeColor','none','FaceColor','b','FaceAlpha',fa);
[pp1,xp1,bwp1] = ksdensity(results.theta1(burn_in:end,:));
plot(xp1,pp1,'color','b','linewidth',2);
% Plot true theta1
% set(gca,'YLim',[0,4.5]);
ylims = get(gca,'YLim');
plot([theta1 theta1],ylims,'--','Color',lcol,'LineWidth',1.5);
fillunder(xp1,pp1,'b',falpha);
set(gca,'YLim',ylims);
% Put a legend on it
lg1 = legend('Prior','Posterior','Estimate');
% title('Prior and posterior distributions of \theta_1');
flushLegend(lg1,'northeast');
xlabel('Elastic modulus');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 for KOH
% histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
%     'EdgeColor','b',...
%     'LineWidth',2,...
%     'FaceAlpha',fa,...
%     'BinWidth',1,...
%     'DisplayStyle','stairs');
[pp2,xp2,bwp2] = ksdensity(results.theta2(burn_in:end,:),'Bandwidth',0.4);
plot(xp2,pp2,'color','b','linewidth',2);
fillunder(xp2,pp2,'b',falpha);

% Put a legend on it
lg2 = legend('Prior','Posterior');
% title('Prior and posterior distributions of \theta_2');
xlabel('Gain');
flushLegend(lg2,'northeast');
set(f,'color','white');

% Set suptitle and fix positions 
suptitle('Prior and posterior distributions for dynamic vibration system');
f.Children(2).Position(2) = .2 ; f.Children(2).Position(4) = .625 ; 
f.Children(5).Position(2) = .2 ; f.Children(5).Position(4) = .625 ; 
flushLegend(lg2,'northeast');
axes(f.Children(5)); flushLegend(lg1,'northwest');
% TODO fix title overlap scale

% Save it:
savestr = ...
sprintf(['FIG_DVS_DCTO_input_posteriors']);
set(f,'PaperPositionMode','auto')
% print(f,savestr,'-depsc','-r600')

%% Get prior and post. predictive distributions
clc ; clearvars -except dpath ; close all ;

% Load previously gathered results
loadloc = [dpath,'dual_calib\DVS_application\data\',...
        '\2020-09-30_dvs_dcto_results'];
load(loadloc,'results');
burn_in=results.settings.burn_in; 

% Let user know we're working here
fprintf('\nBeginning estimation of prior and posterior outputs\n');

% Define inputs mins and ranges 
theta1 = 6.2e10;
xmin = results.settings.min_x;
xrange = results.settings.range_x;
t1min = results.settings.min_t1;
t1range = results.settings.range_t1;
t2min = results.settings.min_t2;
t2range = results.settings.range_t2;
ystd = results.settings.std_y;
ymean = results.settings.mean_y;
burn_in = results.settings.burn_in;

% Get theta1, theta2, rho, lambda, observations and targets from dcto
t1 = results.theta1(burn_in+1:end);
t1_01 = (t1-t1min)./t1range;
t2 = results.theta2(burn_in+1:end);
t2_01 = (t2-t2min)./t2range;
rho = results.obs_rho(burn_in+1:end,:);
lambda = results.obs_lambda(burn_in+1:end,:);
% Recover observations and target outcomes
obs_x_01  = results.settings.obs_x;
obs_y_std  = results.settings.obs_y;
obs_t2_01 = results.settings.obs_t2;
des_y_std  = results.settings.des_y;
m = size(t1,1);
prior_mean = results.settings.mean_obs;

% Set desired observations
n = 3 ; % Number of points to use for integration
x_01 = linspace(0,1,n)';% Get points 
x = x_01 * xrange + xmin; 

% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-5*eye(size(X)); % adds nugget for computational stablty
sig2=0.05; % known observation variance
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    prior_mean(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    (y - prior_mean(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;

% Get computer model output for each draw from the posterior (at x),
% and also get true output
comp_model_output=emulator_mean(results,repmat(x,m,1),...
    [repelem(t1,n,1),repelem(t2,n,1)]);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
comp_model_output = reshape(comp_model_output,n,m)';

% Add posterior mean of observation discrepancy
discrep_gp_post_means_std = ...
    comp_model_output * 0 ;%Pre-allocate space
discrep_gp_post_sds_std = ...
    comp_model_output * 0 ; % Pre-allocate space
obs_disc_std = zeros(m,size(obs_x_01,1)) ; % Pre-allocate space
for idx = 1:m
    rho_dcto = rho(idx,:) ; 
    lambda_dcto = lambda(idx,:) ; 
    t1_dcto = t1(idx,:); 
    t2_dcto = t2(idx,:);
    % obs_disc_std holds the discrepancy between the observed values of y
    % and the computer model output for each draw of theta_1, on the
    % standardized scale.
    obs_disc_std(idx,:)= ...
        obs_y_std' - ...
        (emulator_mean(results,...
            obs_x_01*xrange+xmin,...
            [t1_dcto*ones(size(obs_t2_01)),obs_t2_01*t2range+t2min])-ymean)./ystd;
    % Gather discrepancy mean on standardized scale:
    d_dcto=updated_mean(obs_disc_std(idx,:)',obs_x_01,obs_t2_01,...
        x_01,(t2_dcto-t2min)./t2range,rho_dcto,lambda_dcto);
    discrep_gp_post_means_std(idx,:) = d_dcto ;
    % Now get standard deviations of d:
    d_std_cov_dcto = ...
        updated_cov(obs_x_01,obs_t2_01,x_01,(t2_dcto-t2min)./t2range,rho_dcto,lambda_dcto) ; 
    d_std_sds_dcto = sqrt(diag(d_std_cov_dcto)+0.0) ; % Could add sig2 here
    discrep_gp_post_sds_std(idx,:) = d_std_sds_dcto ; 
    if mod(idx,1000)==0 ; disp(m-idx) ; end
end

% Add discrepancy means to computer model output, rescale sds
discrep_gp_post_sds = ...
    discrep_gp_post_sds_std * ystd ;% Put on orig scale
discrep_gp_post_means = ...
    discrep_gp_post_means_std * ystd;% Put on orig scl
posterior_preds = comp_model_output + discrep_gp_post_means;
posterior_lb = ...
    posterior_preds - 2*discrep_gp_post_sds; % lower bound
posterior_ub = ...
    posterior_preds + 2*discrep_gp_post_sds; % upper bound

% Generate samples of theta1,theta2 from prior distribution
theta1 = rand(m,1) * t1range + t1min;
theta2 = rand(m,1) * t2range + t2min;

% Get predicted output at each x point
prior_model_output = emulator_mean(results,repmat(x_01,m,1),...
    [repelem(theta1,n,1),repelem(theta2,n,1)]);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
prior_model_output = reshape(prior_model_output,n,m)';

% Help function
fillunder = @(x,y,color,falpha) ...
    fill([x(1) x x(end) fliplr(x) 0],...
        [0 y 0 0*y 0],color,'EdgeColor','none','FaceAlpha',falpha);

% Show posterior predictive distribution of average output,
% and include prior predictive distribution
f=figure('pos',[10 10 800 215]);
falpha=.1;
for ii = 1:n
    posterior_distro = ...
        mean(normpdf(linspace(0,1),...
        posterior_preds(:,ii),discrep_gp_post_sds(:,ii)));
    subplot(1,n,ii);
%     histogram(prior_model_output(:,ii),'Normalization','pdf',...
%         'EdgeColor','none','FaceColor','g','FaceAlpha',.5); hold on;
%     histogram(prior_model_output(:,ii),'Normalization','pdf',...
%         'DisplayStyle','stairs','EdgeColor','k','HandleVisibility','off');
    [ypp,xpp,bwp]=ksdensity(prior_model_output(:,ii));
    plot(xpp,ypp,...
        'color','g',...
        'LineWidth',2,...
        'LineStyle',':')
    hold on;
%     area(linspace(0,1),posterior_distro,...
%         'FaceColor','b',...
%         'FaceAlpha',.1,...
%         'EdgeColor','b',...
%         'LineWidth',2); 
    plot(linspace(0,1),posterior_distro,...
        'color','b',...
        'LineWidth',2);
    
    % Plot NSGA-II result
    nsga_results = unique(results.nsga_result.final_obj,'rows');
    xline(nsga_results(ii),'Linewidth',2,'color','r','Linestyle','--');
    
    fillunder(xpp,ypp,'g',falpha);
    fillunder(linspace(0,1),posterior_distro,...
        'b',falpha);
    title(sprintf('mass = %gkg',x(ii)));
    xlim([0,1.2]);
    set(gca,'YTick',[]);
    
%     set(gca,'YLim',ylims);
    if ii==1 
        lg=legend('Prior','Posterior','NSGA-II','Location','northeast'); 
        legend boxoff;
    end
end
% Set title

suptitle('Prior and posterior damping ratio at various oscillator masses');


% Save it:
set(f,'Color','w');
savestr = ...
sprintf(['FIG_DVS_DCTO_prior_and_posterior_output']);
set(f,'PaperPositionMode','auto')
% print(f,savestr,'-depsc','-r600')