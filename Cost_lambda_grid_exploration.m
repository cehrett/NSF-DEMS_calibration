% Workspace for grid investigation of Cost_lambda

clc; clear all; close all;

% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);

%% User defined values
m = 10; % Grid size
M = 6e3; % MCMC length
desired_obs = [ 0 0 ];
which_outputs = [ 1 1 0 ] ; % Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
settings.burn_in=2000;
% trivar_output_settings is useful for getting output at post mean
trivar_output_settings = MCMC_settings(M,[0 0 0],[1 1 1]);
trivar_output_settings.burn_in=2000;
% Changed the MCMC_settings file (new dummy vector), so to use the results
% I recorded earlier, I need to modify the settings file here for backwards
% compatibility:
trivar_output_settings.sim_xt = [ ...
    (2- 2 * trivar_output_settings.sim_xt(:,1) - ...
    trivar_output_settings.sim_xt(:,2) )/2 ...
    trivar_output_settings.sim_xt(:,3:end)];
% Also need old omega,rho,lambda
trivar_output_settings.omega  = [0.655344235568109   0.931941001705886];
trivar_output_settings.rho    = [0.960653924901867   0.991953924787049];
trivar_output_settings.lambda = [0.017385994893994                    ];
% And, need old version of obs_x
trivar_output_settings.obs_x = [ ...
    (2- 2 * trivar_output_settings.obs_x(:,1) - ...
    trivar_output_settings.obs_x(:,2) )/2 ...
    trivar_output_settings.obs_x(:,3:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop through grid of Cost_lambda vals
gridvals = linspace(0,100,m); % Get grid values
results = cell(m,1); % This will store results
for ii = 1:4%1:m
    
    settings.Cost_lambda = gridvals(ii);
    
    % Joint prop for theta, joint prop for obs var, prior on obs var
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    
    % Collect results
    post_mean_out = em_out(samples,trivar_output_settings);
    result = struct('samples',samples,...
        'sigma2',sigma2_rec,...
        'Sigma',Sigma,...
        'desired_obs',desired_obs,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'settings',settings);
    
    results{ii+1} = result;
    
    save([dpath,'stored_data\'...
        'results_Cost_lambda_grid_exploration'],...
        'results');
    
end

save([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
    'results');

%% Figures
load([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
    'results');

% Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
% well as upper and lower .05 quantiles
m=size(results,1); % Store number of target cost_lambdas
cost_lambda = zeros(m,1); % This will store cost_lambdas
alpha = 0.05; % Set quantile
pmo = zeros(m,3); % This will store posterior mean output of emulator
pdo = zeros(m,3); % ``'' median output
pso = zeros(m,3); % ``'' 2 standard deviations
plo = zeros(m,3); % ``'' lower (alpha) quantile
puo = zeros(m,3); % ``'' upper (alpha) quantile
for ii = 1:m % This loop populates the above arrays
    pmo(ii,:) = results{ii}.post_mean_out;
    pdo(ii,:) = quantile(results{ii}.model_output.by_sample,0.5);
    pso(ii,:) = 2 * mean(results{ii}.model_output.sds);
    plo(ii,:) = quantile(results{ii}.model_output.by_sample,alpha);
    puo(ii,:) = quantile(results{ii}.model_output.by_sample,1-alpha);
    cost_lambda(ii) = results{ii}.Cost_lambda;
end
% Now we break the arrays up each into 3 vectors, one for each output
post_cost_mean = pmo(:,3);
post_defl_mean = pmo(:,1);
post_rotn_mean = pmo(:,2);
post_cost_median = pdo(:,3);
post_defl_median = pdo(:,1);
post_rotn_median = pdo(:,2);
post_cost_sd = pso(:,3);
post_defl_sd = pso(:,1);
post_rotn_sd = pso(:,2);
post_cost_lq = plo(:,3);
post_cost_uq = puo(:,3);
post_defl_lq = plo(:,1);
post_defl_uq = puo(:,1);
post_rotn_lq = plo(:,2);
post_rotn_uq = puo(:,2);

% Begin the figure
h=figure('rend','painters','pos',[10 10 1200 400]);
x = 0:1:100;
subplot(1,3,1)
% Get main curve
pcost = pchip(cost_lambda,post_cost_median,x);
% Get upper and lower 0.05 quantiles curves
pcostuq = pchip(cost_lambda,post_cost_uq,x);
pcostlq = pchip(cost_lambda,post_cost_lq,x);
f=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost_lambda,post_cost_median,'or',...
    x,pcost,'-r',...
    x,pcostuq,'-k',...
    x,pcostlq,'-k');
% Plot 2sd errbar
% errorbar(cost_lambda,post_cost_mean,post_cost_sd,'ob'); 
xl1=xlabel('\lambda_c_o_s_t');

ylabel('Cost');

subplot(1,3,2)
% Get main curve
pdefl = pchip(cost_lambda,post_defl_median,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(cost_lambda,post_defl_uq,x);
pdefllq = pchip(cost_lambda,post_defl_lq,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost_lambda,post_defl_median,'or',...
    x,pdefl,'-r',...
    x,pdefluq,'-k',...
    x,pdefllq,'-k');
xl2=xlabel('\lambda_c_o_s_t');
ylabel('Deflection');

subplot(1,3,3)
% Get main curve
protn = pchip(cost_lambda,post_rotn_median,x);
% Get upper and lower 0.05 quantiles curves
protnuq = pchip(cost_lambda,post_rotn_uq,x);
protnlq = pchip(cost_lambda,post_rotn_lq,x);
f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost_lambda,post_rotn_median,'or',...
    x,protn,'-r',...
    x,protnuq,'-k',...
    x,protnlq,'-k');
xl3=xlabel('\lambda_c_o_s_t');
ylabel('Rotation');

% Put a main title over anything, and fix any misplaced labels etc
suptitle(['Performance metrics vs. \lambda_c_o_s_t,',...
    ' with upper/lower ',num2str(alpha),' quantile bands']); 
p = get(xl1,'position');
set(xl1,'position',p + [0 6 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.003 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.00045 0])

%saveas(h,'FIG_cost_lambda.png');

%% Get confidence intervals at each Cost_lambda
m=size(results,1); % Get total number of distinct Cost_lambdas
n = 0; % Set sample size for estimating mean, sd at each point;
       % n=0 means use all samples (computationally expensive!)

% These will store the results, sd's and means at each Cost_lambda
outputs = cell(m,1);
intervals = zeros(m,3);
means = zeros(m,3);

% Loop: get the emulator output at each sample drawn in each MCMC
for ii = 2:m
    fprintf('Step %d/%d\n',ii,m); % Let us know what step we're on
    
    % Get the outputs for the ii^th MCMC chain
    emout = em_out_many(results{ii}.samples,trivar_output_settings,n);
    
    % Record them (in a couple places, for convenience)
    outputs{ii} = emout;
    model_output.by_sample = emout.output_means;
    means(ii,:) = mean(emout.output_means);
    model_output.means = means(ii,:);
    % intervals(ii,:) = std(emout);
    output_gp_sds = emout.output_sds;
    model_output.sds = output_gp_sds;
    % model_output.sds = intervals(ii,:);
    results{ii}.model_output = model_output;
end

% Loop: get the quantiles at each MCMC chain (rather than just relying on
% the standard deviations

save([dpath,'stored_data\'...
     'results_Cost_lambda_grid_exploration'],...
     'results');
