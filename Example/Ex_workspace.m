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
ea = .25       ;  % edge alpha
fa = .75       ;  % face alpha

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

% Add line at posterior mean after below calibration
pmo = mean(samples(settings.burn_in:end,:)) .* settings.input_calib_ranges;
hold on;
plot3([pmo(1) pmo(1)], [pmo(2) pmo(2)], get(gca,'Zlim'), 'k',...
    'LineWidth',6);

%% Get simulation observations
n_cval  = 3 ; % Number of distinct c values to use
n_theta1 = 12 ; % Number of distinct theta1 values to use
n_theta2 = 12 ; % Number of distinct theta2 values to use
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

save([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12'],...
'raw_dat');


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
settings = MCMC_settings (M,desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
    sim_y,which_outputs,Rho_lam_optimum);
settings.doplot = false;
settings.doplot = true;

% MCMC
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);

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

save([dpath,'Example\Ex_results\'...
'2018-05-28_d0_incl_min_cost'],...
'results');

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

save([dpath,'Example\Ex_results\'...
'2018-05-28_d0_incl_min_cost'],...
'results');


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
    
    save([dpath,'Example\Ex_results\'...
        '2018-05-29_cost_grid_results'],...
        'results');

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

% save([dpath,'Example\Ex_results\'...
%     '2018-06-20_discrepancy_full_calib_G50-p25_lambda_prior'],...
%     'results');

