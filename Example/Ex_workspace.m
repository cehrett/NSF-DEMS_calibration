% Simulation workspace

%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);

%% Take a look at surfaces of example function output
clc ; clearvars -except dpath ; close all ; 


theta1=linspace(0,3);
theta2=linspace(0,6);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 , T2 ] = meshgrid(theta1,theta2); 
ths = cat(2,T1',T2');
Theta = reshape(ths,[],2);

% Set other parameters
c = repmat(2,size(Theta,1),1);

% Now get output values
opc = Ex_sim([c Theta]);

% Normalize the outputs
oscl_n = (opc(:,1)-min(opc(:,1))) / range(opc(:,1)) ;
perf_n = (opc(:,2)-min(opc(:,2))) / range(opc(:,2)) ;
cost_n = (opc(:,3)-min(opc(:,3))) / range(opc(:,3)) ;

% Now take a look at the surfaces
oscls = reshape(oscl_n,[],length(theta1));
perfs = reshape(perf_n,[],length(theta1));
costs = reshape(cost_n,[],length(theta1));

ec = 'black' ;  % edge color
ea = .5       ;  % edge alpha
fa = .95       ;  % face alpha

surf(theta2,theta1,oscls,'FaceColor','red','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa);
axis vis3d;

hold on;
surf(theta2,theta1,perfs,'FaceColor','blue','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa); 
axis vis3d;

surf(theta2,theta1,costs,'FaceColor','green','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa); 
axis vis3d;
xlabel('\theta_2'); ylabel('\theta_1'); zlabel('Outcomes');

% Make rotating gif
set(gcf,'color','white')
viewpt = [-27.2667 10.4000];
view(viewpt);
gif('FIG_toy_sim_surfaces.gif','frame',gcf);
nfms = 120;
for ii = 1:nfms
    viewpt = viewpt + [ 360/nfms 0 ];
    view(viewpt);
%     pause(.05)
    gif
end
% saveas(f1,'FIG_nondom_dir_data_vs_est_MCMC_output.png')

% Add line at posterior mean after below calibration
% pmo = mean(samples(settings.burn_in:end,:)) .* settings.input_calib_ranges;
% hold on;
% plot3([pmo(1) pmo(1)], [pmo(2) pmo(2)], get(gca,'Zlim'), 'k',...
%     'LineWidth',6);

%% Get simulation observations
n_cval  = 3 ; % Number of distinct c values to use
n_theta1 = 8 ; % Number of distinct theta1 values to use
n_theta2 = 8 ; % Number of distinct theta2 values to use
cvals  = linspace(1.95,2.05,n_cval)  ; % Get distinct c values
theta1vals = linspace(0,3,n_theta1) ; % Get distinct theta1 values
theta2vals = linspace(0,6,n_theta2) ; % Get distinct theta2 values

% Make a data frame sim_xt = [c theta1 theta2]
[ sc, st1, st2 ] = ndgrid(cvals,theta1vals,theta2vals) ;
sim_xt = [ sc(:) st1(:) st2(:) ] ;

% Get output
sim_y = Ex_sim(sim_xt);

% Package it
raw_dat = struct('sim_xt',sim_xt,'sim_y',sim_y);

% save([dpath,'Example\Ex_results\'...
% '2019-10-24-raw_dat-3-8-8'],...
% 'raw_dat');


%% Run MCMC
clc ; clearvars -except dpath ; close all;
% Load raw data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_xt = raw_dat.sim_xt;
sim_y = raw_dat.sim_y;
clear raw_dat;

% User defined values
M = 2e4;
desired_obs = [0 0 0];
which_outputs = [ 1 1 1 ] ; % Which of oscl, perf, cost
% Calculated optimum for example simulation:
Rho_lam_optimum  = [  0.280981573480363   0.999189406633873...
   0.600440750045477  0.719652153362981   0.102809702497319...
   0.000837772517865 ] ;

% Settings
settings = MCMC_settings (desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
    sim_y,'which_outputs',which_outputs,'Rho_lam_optimum',Rho_lam_optimum);
settings.doplot = false;
settings.doplot = true;

% MCMC
results = MCMC_sigma_prior_joint_prop(settings);

% Get extra info about results and save everything
post_mean_out = em_out(samples,settings)

results = struct('samples',samples,...
    'sigma2',sigma2_rec,...
    'Sigma',Sigma,...
    'init',samples(1,:),...
    'desired_obs',desired_obs,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
    'post_mean_out',post_mean_out,...
    'settings',settings);

% save([dpath,'Example\Ex_results\'...
% '2018-05-28_d0_incl_min_cost'],...
% 'results');

% Add model predictions for each sample point to results
emout = em_out_many(results.samples,results.settings,0);
model_output.by_sample = emout.output_means;
% Then the means
model_output.means = mean(emout.output_means);
model_output.at_post_means = em_out(results.samples,...
    results.settings);
% Then the standard deviations
model_output.sds = emout.output_sds;
% Now package everything up in the results structs
results.model_output = model_output;

% save([dpath,'Example\Ex_results\'...
% '2018-05-28_d0_incl_min_cost'],...
% 'results');


%% Gather results over grid of cost values
clc; clearvars -except dpath; close all;

m=12;
cost_grid = linspace(15,30,m);

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_xt = raw_dat.sim_xt;
sim_y = raw_dat.sim_y;
clear raw_dat;

% User-defined values
M = 1e4+10;
desired_obs = [ 0 0 0 ] ;
which_outputs = [ 1 1 1 ];
Rho_lam_optimum  = [  0.280981573480363   0.999189406633873...
   0.600440750045477  0.719652153362981   0.102809702497319...
   0.000837772517865 ] ;
cost_var = 0.05; % Specifies known cost variance

% Set up struct to catch everything
results = cell(m,1);

for ii = 1:m
    
    % Get new desired_observations
    desired_obs = [ 0 0 cost_grid(ii) ] ;
    
    % Get settings
    settings = MCMC_settings (M,desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
        sim_y,which_outputs,Rho_lam_optimum);
    settings.doplot = true;
    
    % Modify settings to make cost known
    % First, modify the proposal density for obs var.
    settings.proposal.sigma2_prop_density = @(x,s) ...
        [exp(mvnrnd(log([x(1) x(2)]),s(1:2,1:2))) x(3)];
    % Now modify the prior for obs var.
    settings.log_sigma2_prior = @(sigma2) -log(prod(sigma2(1:2)));
    % Now modify the initial obs var.
    settings.sigma2 = [settings.sigma2(1:2) cost_var];
    % Now modify the Metrop.-Hastings correction for drawing obs var.
    settings.log_sig_mh_correction = @(sig_s,sig) ...
        log(prod(sig_s(1:2)))-log(prod(sig(1:2)));
    % We don't want an informative prior on VF, thk., so remove that:
    settings.log_theta_prior = @(theta,Cost_lambda) 0 ;
    % Okay, now tell it we want progress plots during the MCMC
    settings.doplot = true;
    % And set the burn_in to what we want it to be
    settings.burn_in=2000;
    
    % MCMC
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);

    % Get extra info about results and save everything
    post_mean_out = em_out(samples,settings)
    samples_os = samples .* settings.input_calib_ranges + ...
        settings.input_calib_mins;
    result = struct('samples',samples,...
        'samples_os',samples_os,...
        'sigma2',sigma2_rec,...
        'Sigma',Sigma,...
        'desired_obs',desired_obs,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'settings',settings);
    
    results{ii} = result ; 
    
%     save([dpath,'Example\Ex_results\'...
%         '2018-05-29_cost_grid_results'],...
%         'results');

end

load([dpath,'Example\Ex_results\'...
    '2018-05-29_cost_grid_results'],...
    'results');

% Add estimated model outputs for each sample point
% These will store the results, sd's and means at each cost
outputs = cell(m,1);
intervals = zeros(m,3);
means = zeros(m,3);
n=0; 
for ii = 2:12
    fprintf('Step %d/%d\n',ii,m); % Let us know what step we're on
    
%     % Get the outputs for the ii^th MCMC chain
%     emout = em_out_many(results{ii}.samples,results{ii}.settings,n);
%     % Record them (in a couple places, for convenience)
%     % First, all the outputs (one per sample drawn in MCMC)
%     outputs{ii} = emout;
%     model_output.est_by_sample = emout.output_means;
%     % Then the means
%     means(ii,:) = mean(emout.output_means);
%     model_output.est_means = means(ii,:) ;
%     model_output.est_at_post_means = em_out(results{ii}.samples,...
%         results{ii}.settings);
%     % Then the standard deviations
%     output_gp_sds = emout.output_sds;
%     model_output.est_sds = output_gp_sds;
    % Get true output at each sample
    samps = results{ii}.samples_os;
    true_by_sample = Ex_sim([2*ones(size(samps,1),1) samps]) ;
    model_output.true_by_sample=true_by_sample;
    % Now package everything up in the results structs
    results{ii}.model_output = model_output;
    
    save([dpath,'Example\Ex_results\'...
        '2018-05-29_cost_grid_results'],...
        'results');

end

%% Examine the nondominated solutions in cost grid
% Gather all outcomes found in cost grid, using true function to get output
% at sample points
samps = [];
for ii=1:size(results,1)
    ss = results{ii}.samples_os;
    samps = [samps ; Ex_sim([2*ones(size(ss,1),1) ss ]) ] ; 
end

% Find nondominated outcomes
nondoms = nondominated(samps);
% ( Get direct data from Ex_optimization_workspace, allperfs )
%nondom_ap = nondominated(allperfs);
load([dpath,'Example\Ex_results\'...
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');


% Take a look
Circlesize=50;
figure; h1 = scatter3(samps(:,3),samps(:,1),samps(:,2),...
    Circlesize,samps(:,3),'filled','MarkerFaceAlpha',.4);
figure; h1 = scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4);
axis vis3d;
hold on;

% Compare with data obtained directly
scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    Circlesize,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',.2);

%% Examine the solutions in full calibration against true pareto front
% Gather the post burn-in samples:
load([dpath,'Example\Ex_results\'...
'2018-05-17_d0_incl_min_cost'],...
'results');
samps_os = results.samples_os(results.settings.burn_in:end,:);

% Use the emulator to estimate output at these sample points
% (Have to divide the samples into batches, because otherwise the
% covariance matrix is bigger than MATLAB can handle)
% clear samps_os; % We won't need this
% samps = results.samples(results.settings.burn_in:end,:);
% sets = results.settings; sets.burn_in=1;
% clear results; % Free more memory
% samps1 = samps(1:floor(size(samps,1)/2),:);
% samps2 = samps(floor(size(samps,1)/2)+1:end,:);
% clear samps;
% emout1 = em_out_many(samps1,sets,0);
% emout2 = em_out_many(samps2,sets,0);
load([dpath,'Example\Ex_results\'...
    '2018-05-18_full_calib_emulator_output_estimates'],...
    'emout');
% emout.output_means = [ emout1.output_means ; emout2.output_means ] ;
% emout.output_sds   = [ emout1.output_sds   ; emout2.output_sds   ] ;
% save([dpath,'Example\Ex_results\'...
%     '2018-05-18_full_calib_emulator_output_estimates'],...
%     'emout');

% Use the true function to find the output at these sample points
y_samps_true = Ex_sim( [2*ones(size(samps_os,1),1) samps_os]);
% Get nondominated outcomes
nondoms = nondominated(y_samps_true);
% Alternative: use GP estimates rather than true function:
nondoms = nondominated(emout.output_means);
% ( Get direct data from Ex_optimization_workspace, allperfs )
%nondom_ap = nondominated(allperfs);
load([dpath,'Example\Ex_results\'...
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');

% Take a look
Circlesize=50;
%figure; h1 = scatter3(y_samps(:,3),y_samps(:,1),y_samps(:,2),...
%    Circlesize,'b','filled','MarkerFaceAlpha',.4);
figure; h1 = scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4);
axis vis3d;
hold on;

% Compare with data obtained directly
scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    Circlesize,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',.2);

%% Perform calibration with set total observation variance
% Allowing proportion of observation variance allocated to each output to
% vary throughout the chain.
clc ; clearvars -except dpath ; close all;

% Load raw data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_xt = raw_dat.sim_xt;
sim_y = raw_dat.sim_y;
clear raw_dat;

% User defined values
M = 2e4;
desired_obs = [0 0 0];
which_outputs = [ 1 1 1 ] ; % Which of oscl, perf, cost
% Calculated optimum for example simulation:
Rho_lam_optimum  = [  0.280981573480363   0.999189406633873...
   0.600440750045477  0.719652153362981   0.102809702497319...
   0.000837772517865 ] ;

% Settings
% settings = MCMC_settings(M,desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
%     sim_y,which_outputs,Rho_lam_optimum);
settings = MCMC_settings(desired_obs,sim_xt(:,1),sim_xt(:,2:3),sim_y,...
    'M',200,'ObsVar','STOV');
settings.doplot = false;
settings.doplot = true;

% Alter settings to perform calibration with set total obs. variance
settings.sigma2 = 10; % Total observation variance
settings.log_sigma2_prior = @(x) 0 ; % No prior on total variance
settings.init_sigma2_divs = [ 1/3 2/3 ] ; % Initial weights
settings.log_sig_mh_correction = @(x,s) 0 ; % No MH correction for sigma2
settings.proposal.sigma2_prop_density = ...
    @(x,s) x + rand(1,2) * s - s/2;
settings.proposal.Sigma_sig = .1;

% MCMC
[samples,sigma2_rec,Sigma] = MCMC_set_total_obs_var(settings);

% Get extra info about results and save everything
post_mean_out = em_out(samples,settings)
samples_os = samples .* [3 6 ];

results = struct('samples',samples,...
    'sigma2',sigma2_rec,...
    'Sigma',Sigma,...
    'desired_obs',desired_obs,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
    'post_mean_out',post_mean_out,...
    'settings',settings);

% save([dpath,'Example\Ex_results\'...
% '2018-05-29_set_obs_var_d0'],...
% 'results');

% Add model predictions for each sample point to results
emout = em_out_many(results.samples,results.settings,0);
model_output.by_sample_est = emout.output_means;
model_output.by_sample_true = Ex_sim(...
    [2 * ones(size(samples_os,1),1) samples_os ]);
% Then the means
model_output.means_est = mean(emout.output_means);
model_output.at_post_means_est = em_out(results.samples,...
    results.settings);
% Then the standard deviations
model_output.sds = emout.output_sds;
% Now package everything up in the results structs
results.model_output = model_output;

% save([dpath,'Example\Ex_results\'...
% '2018-05-29_set_obs_var_d0'],...
% 'results');

%% Cost grid calibration using set total observation variance
clc; clearvars -except dpath; close all;

m=12;
cost_grid = linspace(15,30,m);

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_xt = raw_dat.sim_xt;
sim_y = raw_dat.sim_y;
clear raw_dat;

% User-defined values
M = 1e4;
desired_obs = [ 0 0 0 ] ;
which_outputs = [ 1 1 1 ];
Rho_lam_optimum  = [  0.280981573480363   0.999189406633873...
   0.600440750045477  0.719652153362981   0.102809702497319...
   0.000837772517865 ] ;
cost_var = 0.05; % Specifies known cost variance

% Set up struct to catch everything
results = cell(m,1);

for ii = 1:m
    
    % Get new desired_observations
    desired_obs = [ 0 0 cost_grid(ii) ] ;
    
    % Get settings
    settings = MCMC_settings (M,desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
        sim_y,which_outputs,Rho_lam_optimum);
    settings.doplot = true;
    
    % Modify settings to make cost known
    % We don't want an informative prior on VF, thk., so remove that:
    settings.log_theta_prior = @(theta,Cost_lambda) 0 ;
    % Okay, now tell it we want progress plots during the MCMC
    settings.doplot = true;
    % And set the burn_in to what we want it to be
    settings.burn_in=2000;
    
    % Alter settings to perform calibration with set total obs. variance
    settings.sigma2 = 50; % Total observation variance
    settings.log_sigma2_prior = @(x) 0 ; % No prior on total variance
    % Initial weights - this gives cost obs var cost_var
    settings.init_sigma2_divs = [ 1/2 1-cost_var/settings.sigma2 ] ; 
    settings.log_sig_mh_correction = @(x,s) 0 ; % No MH corr'n for sigma2
    settings.proposal.sigma2_prop_density = ...
        @(x,s) [ x(1) + rand * s - s/2 , x(2) ] ;
    settings.proposal.Sigma_sig = .1;

    % MCMC
    [samples,sigma2_rec,Sigma] = MCMC_set_total_obs_var(settings);

    % Get extra info about results and save everything
    post_mean_out = em_out(samples,settings)
    samples_os = samples .* settings.input_calib_ranges + ...
        settings.input_calib_mins;
    result = struct('samples',samples,...
        'samples_os',samples_os,...
        'sigma2_weights',sigma2_rec,...
        'Sigma',Sigma,...
        'desired_obs',desired_obs,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'settings',settings);
    
    results{ii} = result ; 
    
%     save([dpath,'Example\Ex_results\'...
%         '2018-05-30_cost_grid_with_set_total_obs_var'],...
%         'results');

end

load([dpath,'Example\Ex_results\'...
    '2018-05-30_cost_grid_with_set_total_obs_var'],...
    'results');

% Add estimated model outputs for each sample point
% These will store the results, sd's and means at each cost
outputs = cell(m,1);
intervals = zeros(m,3);
means = zeros(m,3);
n=0; 
for ii = 2:12
    fprintf('Step %d/%d\n',ii,m); % Let us know what step we're on
    
%     % Get the outputs for the ii^th MCMC chain
%     emout = em_out_many(results{ii}.samples,results{ii}.settings,n);
%     % Record them (in a couple places, for convenience)
%     % First, all the outputs (one per sample drawn in MCMC)
%     outputs{ii} = emout;
%     model_output.est_by_sample = emout.output_means;
%     % Then the means
%     means(ii,:) = mean(emout.output_means);
%     model_output.est_means = means(ii,:) ;
%     model_output.est_at_post_means = em_out(results{ii}.samples,...
%         results{ii}.settings);
%     % Then the standard deviations
%     output_gp_sds = emout.output_sds;
%     model_output.est_sds = output_gp_sds;
    % Get true output at each sample
    samps = results{ii}.samples_os;
    true_by_sample = Ex_sim([2*ones(size(samps,1),1) samps]) ;
    model_output.true_by_sample=true_by_sample;
    % Now package everything up in the results structs
    results{ii}.model_output = model_output;
    
%     save([dpath,'Example\Ex_results\'...
%         '2018-05-30_cost_grid_with_set_total_obs_var'],...
%         'results');

end

%% Calibration using discrepancy function
clc; clearvars -except dpath; close all;

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
desired_obs = [ 0 0 0 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant');

results = MCMC_discrepancy(settings);

% save([dpath,'Example\Ex_results\'...
%     '2018-06-19_discrepancy_full_calib_G25-1_lambda_prior'],...
%     'results');

% save([dpath,'Example\Ex_results\'...
%     '2018-06-19_discrepancy_full_calib_G5-5_lambda_prior'],...
%     'results');

% load([dpath,'Example\Ex_results\'...
%     '2018-06-20_discrepancy_full_calib_G50-p25_lambda_prior'],...
%     'results');

%% STOV calibration using true function (not emulator)
clc; clearvars -except dpath; close all;

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
desired_obs = [ 0 0 0 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'M',2e4,'ObsVar','STOV');

results = MCMC_set_total_obs_var_true_fn(settings);

%save([dpath,'Example\Ex_results\'...
%     '2018-06-27_STOV_true_fn'],...
%     'results');

%% Calibration using discrepancy using true function
clc; clearvars -except dpath; close all;

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
desired_obs = [ 0 0 0 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant');

results = MCMC_discrepancy_true_fn(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);


% save([dpath,'Example\Ex_results\'...
%     '2018-06-28_discrepancy_true_fn_G50-p5_lambda_prior'],...
%     'results');

%% Take a look at predictions from discrepancy calib using true function
clc ; clearvars -except dpath ; close all ;

% load results
load([dpath,'Example\Ex_results\'...
    '2018-06-28_discrepancy_true_fn_G50-p5_lambda_prior'],...
    'results');
burn_in  = results.settings.burn_in;
ngridpts = 21;

% Set prediction points
xpred=((linspace(1.95,2.05,ngridpts)-results.settings.input_cntrl_mins)/...
    results.settings.input_cntrl_ranges)';

pred_pts = [ repmat([1 0],size(xpred,1),1) xpred ; ...
             repmat([0 1],size(xpred,1),1) xpred ; ...
             repmat([0 0],size(xpred,1),1) xpred ] ;

% First check just the first sample
delta_draw = results.delta_samps(burn_in+1,:)   ;
omega = delta_draw(1:3);
rho   = []             ;
lambda= delta_draw(4)  ;
theta_draw = results.samples_os(burn_in+1,:)    ;

% Now get the training points (ie the desired data points)
sim_dat_input = results.settings.obs_x;
% Get the inputs back into original scale
obs_x_os = sim_dat_input(:,3) * results.settings.input_cntrl_ranges + ...
    results.settings.input_cntrl_mins;
obs_x_os = unique(obs_x_os);
% Get the true output at obs_x_os
true_y_os = Ex_sim([obs_x_os repmat(theta_draw,size(obs_x_os,1),1)]);
% Convert this to standardized scale
true_y = (true_y_os - results.settings.output_means')./...
    results.settings.output_sds';
true_y = true_y(:);
% Get "observed" output at obs_x_os
des_obs = results.settings.y;
% Get discrepancy value at observed points
sim_dat_output = des_obs - true_y;

% Get emulator output
em = emulator(sim_dat_input,sim_dat_output,pred_pts,omega,rho,lambda,...
    0,0,true);

% take a look
disc_outs = reshape(em.mu,size(em.mu,1)/3,[]);
sim_dat_output = reshape(sim_dat_output,size(sim_dat_output,1)/3,[]);
sdo_spread = nan(size(xpred,1),size(sim_dat_output,2));
sdo_spread(1,:) = sim_dat_output(1,:);
sdo_spread(floor(ngridpts/2)+1,:) = sim_dat_output(2,:);
sdo_spread(ngridpts,:) = sim_dat_output(3,:);
disp([xpred disc_outs sdo_spread])

%% Get desired observation specified distance from estimated Pareto front
% On standardized scale
clc ; clearvars -except dpath ; close all ;

% Set desired observation (on original scale)
des_obs = [ 0 0 0 ];

% Set desired distance from Pareto front for new des obs found below
spec_dist = 15.5;

%%% Get points estimating the PF
% Use results from set total observation variance method:
load([dpath,'Example\Ex_results\'...
'2018-05-28_d0_incl_min_cost'],...
'results');
samps_os = results.samples_os(results.settings.burn_in+2:end,:);

% load([dpath,'Example\Ex_results\'...
%     '2018-05-18_full_calib_emulator_output_estimates'],...
%     'emout');
emout.output_means = results.model_output.by_sample_est(...
    results.settings.burn_in+2:end,:);

% Use GP estimates rather than true function:
[nondoms,ndidx] = nondominated(emout.output_means);
nondom_inputs = samps_os(ndidx,:);

% Put the model output estimates on the standardized scale
nondoms_std = (nondoms - mean(results.settings.output_means'))./...
    mean(results.settings.output_sds');

% Put des_obs into standardized scale
des_obs_std = (des_obs - mean(results.settings.output_means'))./...
    mean(results.settings.output_sds');

%%% Get distances from desired observation to Pareto front
dists = sqrt(sum((des_obs_std - nondoms_std).^2,2)) ;
[mindist,mindex] = min(dists);
% Get closest point in PF to des obs (closest in std scale, pt in orig sc)
pf_optim = nondoms(mindex,:);
pf_optim_std = nondoms_std(mindex,:);

% Take a look at the point we chose
scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2)); hold on;
line([0 pf_optim(3)],[0 pf_optim(1)],[0 pf_optim(2)],...
    'Color','r');
dists_os=sqrt(sum((des_obs - nondoms).^2,2)) ;
% Also check out closest point on original scale, in green
[mindist_os,mindex_os]=min(dists_os);
line([0 nondoms(mindex_os,3)],[0 nondoms(mindex_os,1)],...
    [0 nondoms(mindex_os,2)],'Color','g');

%%% Find a point that is set distance (eg 1) from pareto front in the
%%% direction of the original desired observation
%X%dirvec = pf_optim - des_obs; dirvec_normd = dirvec / norm(dirvec);
dirvec = pf_optim_std - des_obs_std; 
dirvec_normd = dirvec/norm(dirvec);
%X%des_obs_new = pf_optim - spec_dist * dirvec_normd ; 
des_obs_new = pf_optim_std - spec_dist * dirvec_normd ;
des_obs_new_os = des_obs_new .* ...
    mean(results.settings.output_sds') + ...
    mean(results.settings.output_means')

% Take a look at the new des obs
%X%plot3(des_obs_new(3),des_obs_new(1),des_obs_new(2),'ro');
plot3(des_obs_new_os(3),des_obs_new_os(1),des_obs_new_os(2),'go');

%% Perform calibration with discrep marginal var specified, using true fn
% Also use intelligently chosen desired observation that is the specified
% distance from the estimated pareto front. That des obs is found elsewhere
clc ; clearvars -except dpath ; close all; 
% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
desired_obs = [0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; % 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',false,'M',1e4,'ObsVar','Constant',...
    'Rho_lam_optimum',0); % Setting Rho_lam_optimum=0 makes it get MLEs

% Modify settings to use constant lambda_delta
% settings.lambda_delta_init = 1/64; % 1; % 
settings.log_lambda_delta_prior = @(ld) ld;
settings.proposal.lambda_prop_density = @(x,s) x;
settings.proposal.log_mh_correction_ld = @(ld_s,ld) 0 ;

% results = MCMC_discrepancy_true_fn(settings);
results = MCMC_discrepancy(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

% save([dpath,'Example\Ex_results\'...
%     '2019-10-16_no_discrep_no_PCTO_with_emulator'],...
%     'results');

%% Given desired observation, find true optimum of toy sim problem using dd
clc ; clearvars -except dpath ; close all;

% load true Pareto front
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get standardized version
% First load some calib results, just to get the settings
load([dpath,'Example\Ex_results\'...
    '2019-10-16_no_discrep_with_PCTO_with_emulator'],...
    'results');
cntrl_mins = results.settings.input_cntrl_mins   ;
cntrl_rngs = results.settings.input_calib_ranges ;
calib_mins = results.settings.input_calib_mins   ;
calib_rngs = results.settings.input_calib_ranges ;
output_mns = results.settings.output_means       ;
output_sds = results.settings.output_sds         ;

% Set desired observation
des_obs_os = [ 0 0 0 ] ; % on original scale
des_obs    = (des_obs_os - mean(output_mns,2)')./mean(output_sds,2)';

% Find closest point in Pareto front
true_pf_obs_os = ctheta_output_nondom(:,4:6); % Get outputs of true pf, os
true_pf_obs = (true_pf_obs_os - mean(output_mns,2)')./mean(output_sds,2)';
[m,idx] = min(sum((true_pf_obs - des_obs).^2,2));
optim = ctheta_output_nondom(idx,:);
optim_calib = optim(2:3);

% Take a look against a heatmap
calib_heatmap(des_obs_os);
hold on;
plot(optim_calib(1),optim_calib(2),'.g');
burn_in = results.settings.burn_in;
scatter(results.samples_os(burn_in:end,1),...
    results.samples_os(burn_in:end,2));
plot(mean(results.samples_os(burn_in:end,1)),...
    mean(results.samples_os(burn_in:end,2)),'.m','MarkerSize',15);   
load([dpath,'Example\Ex_results\'...
    '2019-10-16_no_discrep_no_PCTO_with_emulator'],...
    'results');
scatter(results.samples_os(burn_in:end,1),...
    results.samples_os(burn_in:end,2));
plot(optim_calib(1),optim_calib(2),'.g','MarkerSize',15);   


%% Get farthest distance in parameter space from desired observation
clc ; clearvars -except dpath ; close all ;

% The highest point in the output space is [1 1 30]. So standardize that:
high_obs_os = [1 1 30];
% Load some results just to get the settings for standardizing
load([dpath,'Example\Ex_results\'...
    '2018-07-11_discrepancy_true_fn_set_lambda_delta_1'],...
    'results');
meanobs = results.settings.output_means(:,2)';
sdobs   = results.settings.output_sds(:,2)';
high_obs = (high_obs_os - meanobs)./sdobs
% Now get the distance from desired observation
desired_obs_os = [ 0 0 0 ];
desired_obs = (desired_obs_os - meanobs)./sdobs
dist = sqrt(sum( (high_obs - desired_obs).^2 ) ) 

% Now for a different desired observation
desired_obs_os = [ 0.7130    0.7144   17.9220 ];
desired_obs = (desired_obs_os - meanobs)./sdobs
dist = sqrt(sum( (high_obs - desired_obs).^2 ) ) 

%% Do "preliminary" CDO with true fn and discrep, to estimate Pareto front
clc ; clearvars -except dpath ; close all ;

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
des_obs = [0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(des_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant',...
    'DiscMargPrecProp',@(x,s) exp(mvnrnd(log(x),s)),...
    'DiscMargPrecLogMHCorr',@(sig_s,sig)log(prod(sig_s))-log(prod(sig)));

% Change settings to make lambda_delta prior a little more vague
settings.log_lambda_delta_prior = @(ld)log(gampdf(ld,100,.01));

results = MCMC_discrepancy_true_fn(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

% Save results
load([dpath,'Example\Ex_results\'...
    '2018-12-19_preliminary_cdo_truefn_discrep_do0_ldgam10p1'],...
    'results');

%% Using previous MCMC samps, see how many are within 2sd of des obs, optim
clc ; clearvars -except dpath ; close all ; 

load([dpath,'Example\Ex_results\'...
    '2018-07-11_discrepancy_true_fn_set_lambda_delta_1'],...
    'results');

%%% Get distances of each sample draw from the desired observation
burn_in = results.settings.burn_in + 2;
samps = results.samples_os(burn_in:end,:) ; 
outputs_os = results.model_output.by_sample_true(burn_in:end,:) ; 
% Need outputs on standardized scale:
outputs = (outputs_os - mean(results.settings.output_means'))./...
    mean(results.settings.output_sds');
des_obs = (results.settings.desired_obs - ...
    mean(results.settings.output_means'))./...
    mean(results.settings.output_sds');
desdists = sqrt( sum ( ( outputs - des_obs ).^2, 2 ) ) ; 

%%% Now find out how many are close to des obs
desprop = sum ( desdists < 2 ) / size(desdists,1) 

%%% Now find how many are close to optimum.
% First get the optimum.
[m,i] = min (desdists) ;
optim_input = samps(i,:);
optim_output = outputs(i,:);
optdists = sqrt( sum( (outputs - optim_output).^2, 2) );
optprop = sum( optdists < 1.96 ) / size(optdists,1)

%% Perform calibration without discrepancy function, using emulator
% But make it possible to use discrepancy and true function, for comparison
clc ; clearvars -except dpath ; close all; 
% Load data
% load([dpath,'Example\Ex_results\'...
% '2019-10-17-raw_dat-200obs']);
load([dpath,'Example\Ex_results\'...
'2019-10-16-raw_dat-3-6-6']);

sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% omega, rho, lambda MLEs for 3x6x6 grid
Rho_lambda_optimum = [0.0489 0.9746 0.9990 0.4442 0.0082 0.2638] ; 

% Get settings
desired_obs = [0 0 0]; % [0.7130 0.7144 17.9220]; % 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',false,'M',5e3,'ObsVar','Constant',...
    'Rho_lam_optimum',Rho_lambda_optimum); % 

% Modify settings to use constant lambda_delta
% settings.lambda_delta_init = 1/64; % 1; %Comment out this line if no disc
settings.log_lambda_delta_prior = @(ld) ld;
settings.proposal.lambda_prop_density = @(x,s) x;
settings.proposal.log_mh_correction_ld = @(ld_s,ld) 0 ;

% results = MCMC_discrepancy_true_fn(settings);
results = MCMC_discrepancy(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

% save([dpath,'Example\Ex_results\'...
%     '2019-10-16_no_discrep_no_PCTO_with_emulator'],...
%     'results');

