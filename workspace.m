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

%% Perform calibration using discrepancy function
clc ; clearvars -except dpath ; close all ;

% Load data
load([dpath,'stored_data\'...
'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6);
clear raw_dat;

% Get settings
desired_obs = [ 0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant');

results = MCMC_discrepancy(settings);

save([dpath,'Example\Ex_results\'...
    '2018-07-24_discrepancy_d0'],...
    'results');

%% Perform preliminary CDO to estimate PF
clc ; clearvars -except dpath ; close all ;

% Load data
load([dpath,'stored_data\'...
'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6);
clear raw_dat;

% Get settings
desired_obs = [0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e3,'ObsVar','Constant','burn_in',.5,...
    'DiscMargPrecProp',@(x,s) exp(mvnrnd(log(x),s)),...
    'DiscMargPrecLogMHCorr',@(sig_s,sig)log(prod(sig_s))-log(prod(sig)));

% Change settings to make lambda_delta prior more vague
settings.log_lambda_delta_prior = @(ld)log(exppdf(ld,5));

results = MCMC_discrepancy(settings);

% save([dpath,'stored_data\'...
%     '2018-07-25_discrepancy_d0'],...
%     'results');

%% Given a set of results, get predicted model outputs at each sample pt
clc ; clearvars -except dpath ; close all ;

% Specify path to results
rpath = [dpath,'stored_data\2018-07-27_discrepancy_d-elbow_d-p2'];

% Load results
load(rpath)

samps = results.samples;

% Too big to do all in one go; so split the samples up
samps0 = results.samples(1:results.settings.burn_in+1,:);
samps1 = results.samples(results.settings.burn_in+2:8001,:);
samps2 = results.samples(8002:12001,:);
samps3 = results.samples(12002:16001,:);
samps4 = results.samples(16002:end,:);
% set burn_in to 1 since we already took it out
settings = results.settings; settings.burn_in = 1;
emout0 = em_out_many(samps0,settings,0);
emout1 = em_out_many(samps1,settings,0);
emout2 = em_out_many(samps2,settings,0);
emout3 = em_out_many(samps3,settings,0);
emout4 = em_out_many(samps4,settings,0);

model_output.by_sample_est = [ emout0.output_means ; 
    emout1.output_means ; 
    emout2.output_means ; 
    emout3.output_means ; 
    emout4.output_means ] ;
model_output.sds_by_sample_est = [ emout0.output_sds ; 
    emout1.output_sds ; 
    emout2.output_sds ; 
    emout3.output_sds ; 
    emout4.output_sds ] ;
results.model_output = model_output ; 

% Save the results
save(rpath,...
    'results');

%% Given preliminary CDO, estimate PF and update desired observation
clc ; clearvars -except dpath ; close all ;

%%% Load the preliminary CDO
load([dpath,'stored_data\'...
    '2018-07-25_discrepancy_d0'],...
    'results');

%%% Estimate the PF
[PF_os, PFidx] = nondominated(results.model_output.by_sample_est) ; 

%%% Put PF on standardized scale
omeans = mean(results.settings.output_means');
osds   = mean(results.settings.output_sds'  );
PF     = (PF_os - omeans)./ osds             ;

%%% Find closet point to des_obs
%orig_des_obs = results.settings.desired_obs  ;
orig_des_obs = [.74 .089 100 ] ; % This pt chosen to get at observed elbow
des_obs = (orig_des_obs - omeans)./osds      ;
[m,i] = min( sum( ( PF - des_obs ).^2, 2 ) ) ;
PF_optim = PF(i,:)                           ;

%%% Get new desired obs specified distance from PF in same dir as original
spec_dist = .2                               ;
dirvec_nonnormed = PF_optim - des_obs        ;
dirvec = dirvec_nonnormed/norm(dirvec_nonnormed) ;
des_obs_new = PF_optim - spec_dist * dirvec  ;
des_obs_new_os = des_obs_new .* osds + omeans; 

%%% Take a look
scatter3(PF_os(:,1),PF_os(:,2),PF_os(:,3))   ;
hold on                                      ;
scatter3(orig_des_obs(1),orig_des_obs(2),orig_des_obs(3))          ;
line([orig_des_obs(1) PF_os(i,1)], [orig_des_obs(2) PF_os(i,2)], ...
    [orig_des_obs(3) PF_os(i,3)])                                  ;
scatter3(des_obs_new_os(1),des_obs_new_os(2),des_obs_new_os(3),'g');

%%% Save new desired observation ;
% save([dpath,'stored_data\'...
%     '2018-07-26_elbow_des_obs_d-p2'],...
%     'des_obs_new_os');

%% Use des_obs chosen in preliminary CDO to do calibration with discrep
clc ; clearvars -except dpath ; close all ; 

%%% Load raw data and desired observation
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6);
clear raw_dat;
load([dpath,'stored_data\'...
    '2018-07-26_elbow_des_obs_d-p2']);

% Get settings
desired_obs = des_obs_new_os; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant',...
    'LambdaDeltaInit',1/(.2^2));

results = MCMC_discrepancy(settings);

% save([dpath,'stored_data\'...
%     '2018-07-27_discrepancy_d-elbow_d-p2'],...
%     'results');

%% Get mean output under prior and under posterior
clc ; clearvars -except dpath ; close all ;

%%% Load the results
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');

%%% Get the mean under the prior
mupr = mean(results.settings.output_means');

%%% Get the mean under the posterior
mupo = mean(results.model_output.by_sample_est(...
    results.settings.burn_in+2:end,:));

%%% Display
mupr
mupo

%% Get desired observations for cost grid
clc ; clearvars -except dpath ; close all ;

%%% Specify the distance of the updated desired observations from the PF
spec_dist = 0.2;

%%% Load the results for the Pareto front estimate
load([dpath,'stored_data\'...
    '2018-07-25_discrepancy_d0'],...
    'results');

%%% Estimate the PF
[PF_os, PFidx] = nondominated(results.model_output.by_sample_est) ; 

%%% Put PF on standardized scale
omeans = mean(results.settings.output_means');
osds   = mean(results.settings.output_sds'  );
PF_all     = (PF_os - omeans)./ osds             ;
% Cut out the costs
PF = PF_all;

%%% Set des_obs and find closet point to des_obs
cost_grid_pts = linspace(96,352,20);
orig_des_obs = [ zeros(size(cost_grid_pts,2),2) cost_grid_pts' ] ; 
des_obs = (orig_des_obs - omeans)./osds      ;
PF_optim=zeros(size(des_obs));
for ii = 1 : size(PF_optim,1)
    [m,i] = min( ( PF(:,3) - des_obs(ii,3) ).^2 ) ;
    PF_optim(ii,:) = PF(i,:)                           ;
end
PF_optim_os = PF_optim .* osds + omeans ; 
%plot3(PF_optim_os(:,1),PF_optim_os(:,2),PF_optim_os(:,3),'.r',...
%    'MarkerSize',20)

%%% Now adjust the desired obs so they are close to PF
dirvecs_nonnormed = PF_optim - des_obs ; 
dirvec_norms = sqrt( sum( dirvecs_nonnormed.^2, 2) );
dirvecs = dirvecs_nonnormed ./ dirvec_norms;
des_obs_upd = PF_optim - spec_dist * dirvecs ;
des_obs_upd_os = des_obs_upd .* osds + omeans;

%%% Save new desired observation ;
% save([dpath,'stored_data\'...
%     '2018-08-03_cost_grid_des_obs'],...
%     'des_obs_upd_os');

%% Perform cost_grid calibration
% 2018-08-03
clc ; clearvars -except dpath ; close all ;

%%% Load new desired observation ;
load([dpath,'stored_data\'...
    '2018-08-03_cost_grid_des_obs'],...
    'des_obs_upd_os');

%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6);
clear raw_dat;

%%% Set up for loop over cost grid
m = size(des_obs_upd_os,1);
% results = cell(m,1);

%%% Loop over cost grid, performing CDO at each point
for ii = fliplr(1 : 7)
    %%% Announce what's going on
    fprintf(['\n' repmat('#',1,30) '\n STEP %d of %d, Cost $%3.2f\n' ...
        repmat('#',1,30) '\n\n'],ii,m,des_obs_upd_os(ii,3));
    
    %%% Get settings
    des_obs = des_obs_upd_os(ii,:);
    settings = MCMC_settings(des_obs,sim_x,sim_t,sim_y,...
        'Discrepancy',true,'M',8e3,'ObsVar','Constant',...
        'LambdaDeltaInit',1/(.2^2),'burn_in',2e3/8e3);
    
    %%% Perform calibration
    result = MCMC_discrepancy_costgrid(settings);
    
    %%% Save results
    results{ii} = result;
    save([dpath,'stored_data\'...
    '2018-08-03_cost_grid_discrepancy_results'],...
    'results');

end

% load([dpath,'stored_data\'...
%     '2018-08-03_cost_grid_discrepancy_results'],...
%     'results');

%% Get model output estimates for each point in cost_grid analysis
clc ; clearvars -except dpath ; close all ;

%%% Load the cost_grid results
load([dpath,'stored_data\'...
    '2018-08-03_cost_grid_discrepancy_results'],...
    'results');
m = size(results,1);

%%% Loop through and get the model output estimates
for ii = 1:m 
    
    fprintf(['\n\n' ...
        repmat('#',1,30) '\nSTEP %d\n' repmat('#',1,30) '\n\n'],ii);
    
    res = results{ii};
    samps=res.samples;
    settings=res.settings;
    settings.burn_in = 1 ;
    
    %%% Split the samples up and get output for each subset
    model_output.by_sample_est = [] ; 
    model_output.by_sample_sds = [] ; 
    for jj = 1:4
        subsamps = samps((jj-1)*2000+2-1*(jj==1):jj*2000+1,:);
        subemout = em_out_many(subsamps,settings,0);
        model_output.by_sample_est = [model_output.by_sample_est ; ...
            subemout.output_means ] ;
        model_output.by_sample_sds = [model_output.by_sample_sds ; ...
            subemout.output_sds   ] ;
    end
    
    %%% Save the output estimates to the results
    results{ii}.model_output = model_output;
    results{ii}.post_mean_out = ...
        mean(results{ii}.model_output.by_sample_est(2002:end,:));
    
end

%%% Check for negative variances
minvars = [];
for ii = 1 : m
    minvars = [minvars ; min(min(results{ii}.model_output.by_sample_sds))];
end
min(minvars)

%%% Save the results
% save([dpath,'stored_data\'...
%     '2018-08-03_cost_grid_discrepancy_results'],...
%     'results');

