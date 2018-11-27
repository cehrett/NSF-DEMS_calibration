%%% State-aware new system tinkering workspace

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

%% Define new univariate system
clc ; clearvars -except dpath ; close all ;

% Define control input, calib input and objective function
x = linspace(2,3) ;
t = linspace(0,3) ;
f1 = @(x,t) (t.^(x-1) .* exp(-0.75 * t) + 1) .^ (-1) ;

% Get results on grid 
[X,T] = meshgrid(x,t) ; 
Y = reshape(f1(X(:),T(:)),length(x),length(t)) ;

% Take a look
surf(X,T,Y);

% Get prior mean, sd
output_mean = mean(Y(:));
output_sd   = std(Y(:));
disp(output_mean);
disp(output_sd);

% Looks, as expected, to be the case that theta is optimal at 4/3(x-1).

%% Perform state-aware calibration on univariate system
clc ; clearvars -except dpath ; close all ;

% Set inputs for calibration
desired_obs = 0 ; 
output_mean = 0.6922;
output_sd   = 0.1163;
lambda_delta = 1/(.5)^2;
cntrl_input = linspace(0,1,8);
% f1 = @(x,t) (t.^(x-1) .* exp(-0.75 * t) + 1) .^ (-1) ;
% f2 = @(x,t,varargin) (t.^(x-1) .* exp(-0.75 * t) + 1) .^ (-1) ;
eta = @Univar_objective ; 

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',2,'input_cntrl_range',1,'which_sa',1,...
    'dim_theta1',1,'dim_theta2',0,...
    'input_calib_min',0,'input_calib_range',4,...
    'eta',eta,...
    'output_mean',output_mean,'output_sd',output_sd,...
    'nu_theta_prior_param',100,...
    'nu_delta_prior_param',100,...
    'lambda_theta_hypers',[1, 1],...
    'M',4e4,...
    'burn_in',1/4);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Compare results to truth
res=mean(results.theta1(10001:end,:)) * results.settings.input_calib_range;
truth=4/3 * (linspace(2,3,8) - 1);
figure();
plot(linspace(2,3,8),res,'bo-'); hold on ; plot(linspace(2,3,8),truth) ; 

% Save results
% save([dpath,'state-aware\state-aware_results\'...
%     '2018-11-16_sa_univariate_system'],...
%     'results');

%% Examine marginal posteriors of theta1 at grid of control settings
burn_in = results.settings.burn_in;
cal_rng = results.settings.input_calib_range;
cal_min = results.settings.input_calib_min;

samps = results.theta1(burn_in:end,:) * cal_rng + cal_min;
meantheta = mean(samps);
truetheta = 4/3 * (linspace(2,3,8)-1);
for ii=1:size(samps,2)
    subplot(2,4,ii);
    ksdensity(samps(:,ii));
    xlim([0.5,4]);ylim([0,1.5]);
    hold on;
    plot([meantheta(ii) meantheta(ii)],get(gca,'YLim'));
    plot([truetheta(ii) truetheta(ii)],get(gca,'YLim'));
end