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
    'Discrepancy',true,'M',2e4,'ObsVar','Constant');

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


