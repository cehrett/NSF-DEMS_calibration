%%% State-aware tinkering workspace

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
addpath([dpath,'state-aware']);

%% Perform calibration
clc ; clearvars -except dpath ; close all ;

% Define inputs for calibration:
desired_obs = [0.7130 0.7144 17.9220]; % Chosen to lie dist 1 from est PF
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',1,'input_cntrl_range',2,'which_sa',[1 0]);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Take a look
mean(results.theta1(results.settings.burn_in+1:end,:))
mean(exp(-exp(results.xi(results.settings.burn_in+1:end))))
mean(exp(-exp(results.nu_theta(results.settings.burn_in+1:end,:))))
mean(results.lambda_theta(results.settings.burn_in+1:end))
mean(exp(-exp(results.nu_delta(results.settings.burn_in+1:end,:))))

% Save results
% load([dpath,'state-aware\state-aware_results\'...
%     '2018-10-30_sa_true_fn_ld1_t1calib'],...
%     'results');

%% Perform calibration with the second input state-aware, cntrl rng [1,3]
clc ; clearvars -except dpath ; close all ;

% Define inputs for calibration:
desired_obs = [0.7130 0.7144 17.9220]; % Chosen to lie dist 1 from est PF
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',1,'input_cntrl_range',2,'which_sa',[0 1]);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Take a look
mean(results.theta1(results.settings.burn_in+1:end,:))
mean(exp(-exp(results.xi(results.settings.burn_in+1:end))))
mean(exp(-exp(results.nu_theta(results.settings.burn_in+1:end,:))))
mean(results.lambda_theta(results.settings.burn_in+1:end))
mean(exp(-exp(results.nu_delta(results.settings.burn_in+1:end,:))))

% Save results
% load([dpath,'state-aware\state-aware_results\'...
%     '2018-10-30_sa_true_fn_ld1_t2calib'],...
%     'results');


%% Gather non-state-aware calibration across same control input range [1,3]
clc ; clearvars -except dpath ; close all ; 

% Set inputs for calibration
desired_obs = [ 0.8182    0.6920   17.9234 ] ; 
% desired_obs = [0.7130 0.7144 17.9220];
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',1,'input_cntrl_range',2,...
    'dim_theta1',0,'dim_theta2',2,'which_sa',[0 0]);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn(settings);

% Save results
% load([dpath,'state-aware\state-aware_results\'...
%     '2018-10-29_non-sa_true_fn_ld1'],...
%     'results');

%% Get predictive distribution from existing calibration results
clc ; clearvars -except dpath ; close all ;

% Load relevant calibration results
results_path = [dpath,'state-aware\state-aware_results\'...
    '2018-11-08_non-sa_true_fn_ld1_x23'] ;
load(results_path,'results');

% Extract relevant inputs
xx = results.settings.cntrl_input_with_dum_vars ;
theta1 = results.theta1 ;
theta2 = exp(-exp(results.xi)) ;
input_cntrl_min   = results.settings.input_cntrl_min ; 
input_cntrl_range = results.settings.input_cntrl_range ; 
input_calib_min   = results.settings.input_calib_min ; 
input_calib_range = results.settings.input_calib_range ; 
output_mean       = results.settings.output_mean ;
output_sd         = results.settings.output_sd ;
which_sa          = results.settings.which_sa ;
M                 = results.settings.number_of_iterations ;
dim_output        = results.settings.dim_output ;
dim_cntrl_input   = numel(results.settings.cntrl_input) ;

% Pre-allocate array for predictions
y = nan(M,dim_cntrl_input,dim_output) ;

% Loop through samples
for ii = 1:M
    
    y(ii,:,:) = reshape(Ex_sim_compwise(xx,theta1(ii,:,:),theta2(ii,:),...
        input_cntrl_min,input_cntrl_range,...
        input_calib_min,input_calib_range,...
        output_mean,output_sd,which_sa),...
        dim_cntrl_input,dim_output) ; 
    
    if mod(ii,1000) == 0 disp(ii); end
    
end

% Save results
results.y = y ;
% save(results_path,'results');

%% Get prior output mean and sd under given range of control inputs
clc ; clearvars -except dpath ; close all ;

% Set range of control inputs
input_cntrl_min   = 2 ; 
input_cntrl_range = 1 ;

% Settings for simulation
nsim = 1e5 ; 
cntrl_grid_len = 8 ; 

% Get array of inputs for simulation
theta1 = rand(nsim,1) * 3 ; 
theta2 = rand(nsim,1) * 6 ;
x      = linspace(input_cntrl_min,input_cntrl_min+input_cntrl_range,...
    cntrl_grid_len) ;
x      = repmat(x',nsim,1) ;
inputs = [ x repelem(theta1,cntrl_grid_len,1) ...
    repelem(theta2,cntrl_grid_len,1) ] ; 

% Perform simulation
outputs = Ex_sim(inputs) ; 

% Get mean and sd
output_mean = mean(outputs) 
output_sd   = std(outputs) 

% Also get mean distance from desired observation on standardized scale
%desired_obs = [0 0 0] ;
desired_obs = [0.7587 0.6479 19.1066] ; %[0.7130 0.7144 17.9220] ;
output_mean = [0.9291    0.7182   21.0030] ; %[0.9264 0.7926 21.0015] ;
output_sd   = [0.0550    0.0913    3.2004] ; %[0.0602 0.0900 3.1971] ; 
desired_obs = (desired_obs - output_mean) ./ output_sd ; 
outputs_std = (outputs - output_mean) ./ output_sd ; 
dists = sqrt(sum( (outputs_std - desired_obs).^2, 2 )) ; 
meandist = sum(dists) / length(dists);

% Show them
fprintf('Mean output: ') ; 
fprintf(' %f ',output_mean) ; 
fprintf('\nSd of output: ');
fprintf(' %f ',output_sd);
fprintf('\nMean distance from des_obs on std scale: %f\n',meandist);

%% Compare state-aware (t1 and t2) performance to non-sa performance
clc ; clearvars -except dpath ; close all ; 

% Load non-state-aware calibration results and desired obs
load([dpath,'state-aware\state-aware_results\'...
    '2018-11-08_non-sa_true_fn_ld1_x23'],...
    'results');
results_nsa = results;
desired_obs = results.settings.desired_obs;

% Load state-aware calibration results, t1 as state-aware
load([dpath,'state-aware\state-aware_results\'...
    '2018-12-08_sa_true_fn_ld1_t1calib_x23'],...
    'results');
results_sat1 = results;
clear results;

% Load state-aware calibration results, t2 as state-aware
load([dpath,'state-aware\state-aware_results\'...
    '2018-11-08_sa_true_fn_ld1_t2calib_x23'],...
    'results');
results_sat2 = results;
clear results;

% Convert desired_obs to standardized scale
desired_obs = (desired_obs - results_nsa.settings.output_mean) ./ ...
    results_nsa.settings.output_sd ; 

% Get MSPE of nsa results
diffs_nsa = nan(size(results_nsa.y)) ; 
for ii = 1:3
    diffs_nsa(:,:,ii) = results_nsa.y(:,:,ii) - desired_obs(ii) ; 
end
dists_nsa = sqrt(sum( diffs_nsa.^2, 3)) ;
MSPE_nsa = sum(dists_nsa(:)) / numel(dists_nsa) ; 

% Get MSPE of sat1 results
diffs_sat1 = nan(size(results_sat1.y)) ; 
for ii = 1:3
    diffs_sat1(:,:,ii) = results_sat1.y(:,:,ii) - desired_obs(ii) ; 
end
dists_sat1 = sqrt(sum( diffs_sat1.^2, 3)) ;
MSPE_sat1 = sum(dists_sat1(:)) / numel(dists_sat1) ; 

% Get MSPE of sat1 results
diffs_sat2 = nan(size(results_sat2.y)) ; 
for ii = 1:3
    diffs_sat2(:,:,ii) = results_sat2.y(:,:,ii) - desired_obs(ii) ; 
end
dists_sat2 = sqrt(sum( diffs_sat2.^2, 3)) ;
MSPE_sat2 = sum(dists_sat2(:)) / numel(dists_sat2) ; 

% Compare
fprintf('MSPE for non-state-aware: %f\n',MSPE_nsa); 
fprintf('MSPE for state-aware t1 : %f\n',MSPE_sat1); 
fprintf('MSPE for state-aware t2 : %f\n',MSPE_sat2); 

%% Plot posterior distribution from non-state-aware calibration
clc ; clearvars -except dpath ; close all ; 

% Load results
load([dpath,'state-aware\state-aware_results\'...
    '2018-10-29_non-sa_true_fn_ld1'],...
    'results');
burn_in = results.settings.burn_in ; 
theta1 = exp(-exp(results.xi(burn_in:end,1))) *3; 
theta2 = exp(-exp(results.xi(burn_in:end,2))) *6; 

% Plot
plot( theta1, theta2, 'ko') ;

%% Plot posterior distribution from state-aware theta2 calibration
clc ; clearvars -except dpath ; close all ;

% Load results

burn_in = results.settings.burn_in ; 
theta1 = exp(-exp(results.xi(burn_in:end,1))) *3; 
theta2 = results.theta1(:) * 6;

%% Get heatmap of parameter space, averaged across control input
clc ; clearvars -except dpath ; close all ;

% Set number of points for each calib and control input
M = 500 ;
m = 40 ;

% Set desired observation
desired_obs = [ 0 0 0 ] ;

% Get grid on parameter space and control input space
theta1 = linspace(0,3,M) ; 
theta2 = linspace(0,6,M) ;
x = linspace(2,3,m) ; 
[c,t1,t2] = meshgrid(x,theta1,theta2) ; 
all_inputs = [ c(:) t1(:) t2(:) ] ;
all_outputs = Ex_sim(all_inputs) ; 

% Now get the mean across control inputs
[thetas, ~, Ithetas] = unique(all_inputs(:,2:3),'rows','stable') ;
meanouts = nan(size(thetas,1),3) ; 
for ii = 1:3
    meanouts(:,ii) = accumarray(Ithetas,all_outputs(:,ii),[],@mean);
end

% Get mean, std of meanouts and use them to standardize meanouts & des_obs
output_mean = mean(meanouts) ; 
output_sd   = std(meanouts) ;
meanouts_std = (meanouts - output_mean) ./ output_sd ;
desired_obs_std = (desired_obs - output_mean) ./ output_sd ;

% Now convert meanouts_std to distances from the desired observation
dists = sqrt(sum( (meanouts_std - desired_obs_std).^2, 2)) ;

% Make heatmap
colormap(flipud(jet));
scatter(thetas(:,1),thetas(:,2),2,dists);

% Put true optimum
[~,mindx] = min(dists);
hold on ; 
optim = thetas(mindx,:);
p=plot(optim(1),optim(2),'ok','MarkerSize',7,'MarkerFaceColor','m',...
    'LineWidth',2);

% Add colorbar
colorbar('East');

%% Get updated desired observation close to true optimum
% Note: using brute force here, not preliminary CDO.
clc ; clearvars -except dpath ; close all ;

% Specify distance of new observation from Pareto front
spec_dist = 1;

% Set number of points for each calib and control input
M = 500 ;
m = 40 ;

% Set desired observation
desired_obs = [ 0 0 0 ] ;

% Get grid on parameter space and control input space
theta1 = linspace(0,3,M) ; 
theta2 = linspace(0,6,M) ;
x = linspace(2,3,m) ; 
[c,t1,t2] = meshgrid(x,theta1,theta2) ; 
all_inputs = [ c(:) t1(:) t2(:) ] ;
all_outputs = Ex_sim(all_inputs) ; 

% Now get the mean across control inputs
[thetas, Iin, Ithetas] = unique(all_inputs(:,2:3),'rows','stable') ;
meanouts = nan(size(thetas,1),3) ; 
for ii = 1:3
    meanouts(:,ii) = accumarray(Ithetas,all_outputs(:,ii),[],@mean);
end

% Get mean, std of meanouts and use them to standardize meanouts & des_obs
output_mean = mean(meanouts) 
output_sd   = std(meanouts) 
meanouts_std = (meanouts - output_mean) ./ output_sd ;
desired_obs_std = (desired_obs - output_mean) ./ output_sd ;

% Now convert meanouts_std to distances from the desired observation
dists = sqrt(sum( (meanouts_std - desired_obs_std).^2, 2)) ;

% Get true optimal output
[~,mindx] = min(dists);
optim = meanouts_std(mindx,:);

% Get new desired observation
dirvec = optim - desired_obs_std; 
dirvec_normd = dirvec/norm(dirvec);
des_obs_new = optim - spec_dist * dirvec_normd ;
des_obs_new_os = des_obs_new .* output_sd + output_mean

%% Gather non-state-aware calibration across same control input range [2,3]
clc ; clearvars -except dpath ; close all ; 

% Set inputs for calibration
desired_obs = [ 0.7587    0.6479   19.1066 ] ; 
output_mean = [ 0.9291    0.7182   21.0030 ] ;
output_sd   = [ 0.0550    0.0913    3.2004 ] ;
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',2,'input_cntrl_range',1,...
    'dim_theta1',0,'dim_theta2',2,'which_sa',[0 0],...
    'output_mean',output_mean,'output_sd',output_sd);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn(settings);

% Save results
% load([dpath,'state-aware\state-aware_results\'...
%     '2018-11-08_non-sa_true_fn_ld1_x23'],...
%     'results');

%% Perform calibration with the second input state-aware, cntrl rng [2,3]
clc ; clearvars -except dpath ; close all ;

% Set inputs for calibration
desired_obs = [ 0.7587    0.6479   19.1066 ] ; % ~ dist 1 from est PF
output_mean = [ 0.9291    0.7182   21.0030 ] ;
output_sd   = [ 0.0550    0.0913    3.2004 ] ; 
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',2,'input_cntrl_range',1,'which_sa',[0 1],...
    'dim_theta1',1,'dim_theta2',1,...
    'output_mean',output_mean,'output_sd',output_sd);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Take a look
mean(results.theta1(results.settings.burn_in+1:end,:))
mean(exp(-exp(results.xi(results.settings.burn_in+1:end))))
mean(exp(-exp(results.nu_theta(results.settings.burn_in+1:end,:))))
mean(results.lambda_theta(results.settings.burn_in+1:end))
mean(exp(-exp(results.nu_delta(results.settings.burn_in+1:end,:))))

% Save results
% load([dpath,'state-aware\state-aware_results\'...
%     '2018-11-08_sa_true_fn_ld1_t2calib_x23'],...
%     'results');

%% Perform calibration with the first input state-aware, cntrl rng [2,3]
clc ; clearvars -except dpath ; close all ;

% Set inputs for calibration
desired_obs = [ 0.7587    0.6479   19.1066 ] ; % ~ dist 1 from est PF
output_mean = [ 0.9291    0.7182   21.0030 ] ;
output_sd   = [ 0.0550    0.0913    3.2004 ] ; 
lambda_delta = 1; % Est distance of desired_obs from Pareto front
cntrl_input = linspace(0,1,8); % Control input is on normalized scale

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',2,'input_cntrl_range',1,'which_sa',[1 0],...
    'dim_theta1',1,'dim_theta2',1,...
    'output_mean',output_mean,'output_sd',output_sd);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Take a look
mean(results.theta1(results.settings.burn_in+1:end,:))
mean(exp(-exp(results.xi(results.settings.burn_in+1:end))))
mean(exp(-exp(results.nu_theta(results.settings.burn_in+1:end,:))))
mean(results.lambda_theta(results.settings.burn_in+1:end))
mean(exp(-exp(results.nu_delta(results.settings.burn_in+1:end,:))))

% Save results
save([dpath,'state-aware\state-aware_results\'...
    '2018-11-08_sa_true_fn_ld1_t1calib_x23'],...
    'results');