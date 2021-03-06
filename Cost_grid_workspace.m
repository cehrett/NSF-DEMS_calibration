% Cost grid pareto band workspace
% This is a workspace in which to make a grid of target costs which will be
% specified with extremely low observation variance, to force the MCMC to
% stay at or very near those cost targets. The performance of the samples
% at each grid point will then be plotted along with error bars for 2 sd.

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
desired_def_rot = [ 0 0 ]; % desired deflection and rotation
which_outputs = [ 1 1 1 ] ; % Which of defl, rot, cost
cost_var = 0.05; % Sets constant obs var for cost

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop through grid of Costs
gridvals = linspace(96,350,m); % Get grid values
results = cell(m,1); % This will store results
for ii = 1:m
    
    % Let us know where in the loop we are
    fprintf('Starting run %d/%d:\n\n',ii,m);
    
    % Set desired_obs for this loop
    desired_obs = [desired_def_rot gridvals(ii)];
    
    % Get settings for MCMC
    settings = MCMC_settings(M,desired_obs,which_outputs);
    
    % Modify settings
    % Many settings must be changed to accommodate the fact that we are
    % treating the observation variance for cost as known.
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
    
    % Run the MCMC
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    
    % Collect results
    post_mean_out = em_out(samples,settings);
    result = struct('samples',samples,...
        'sigma2',sigma2_rec,...
        'Sigma',Sigma,...
        'desired_obs',desired_obs,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'settings',settings);
    
    results{ii} = result;
    
    save([dpath,'stored_data\'...
        'results_Cost_grid_exploration'],...
        'results');
    
end

save([dpath,'stored_data\'...
    'results_Cost_grid_exploration'],...
    'results');

load([dpath,'stored_data\'...
    'results_Cost_grid_exploration'],...
    'results');

%% Get output and confidence intervals at each cost
m=size(results,1); % Get total number of distinct costs
n = 0; % Set sample size for estimating mean, sd at each point;
       % n=0 means use all samples (computationally expensive!)

% These will store the results, sd's and means at each cost
outputs = cell(m,1);
intervals = zeros(m,3);
means = zeros(m,3);

% Loop: get the emulator output at each sample drawn in each MCMC
for ii = 1:m
    fprintf('Step %d/%d\n',ii,m); % Let us know what step we're on
    
    % Get the outputs for the ii^th MCMC chain
    emout = em_out_many(results{ii}.samples,results{ii}.settings,n);
    
    % Record them (in a couple places, for convenience)
    % First, all the outputs (one per sample drawn in MCMC)
    outputs{ii} = emout;
    model_output.by_sample = emout.output_means;
    % Then the means
    means(ii,:) = mean(emout.output_means);
    model_output.means = means(ii,:);
    % Then the standard deviations
    output_gp_sds = emout.output_sds;
    model_output.sds = output_gp_sds;
    % Now package everything up in the results structs
    results{ii}.model_output = model_output;
    
    save([dpath,'stored_data\'...
    'results_Cost_grid_exploration'],...
    'results');
    
end

%% Figures
% Collect costs, and posterior mean and sds for defl, rot, as
% well as upper and lower .05 quantiles
m=size(results,1);
cost = zeros(m,1);
cred_level = 90; % Set desired level for credible bands (in %)
alpha = (100-cred_level)/100; % Convert cred_level to alpha level
pmo = zeros(m,3); % This will store posterior mean output of emulator
pdo = zeros(m,3); % ``'' median output
pso = zeros(m,3); % ``'' appropriate multiple of standard deviations
plo = zeros(m,3); % ``'' lower (alpha/2) quantile
puo = zeros(m,3); % ``'' upper (alpha/2) quantile
for ii = 1:m % This loop populates the above arrays
    pmo(ii,:) = results{ii}.post_mean_out;
    pdo(ii,:) = quantile(results{ii}.model_output.by_sample,0.5);
    pso(ii,:) = norminv(1-alpha/2) * mean(results{ii}.model_output.sds);
    plo(ii,:) = quantile(results{ii}.model_output.by_sample,alpha/2);
    puo(ii,:) = quantile(results{ii}.model_output.by_sample,1-alpha/2);
    cost(ii) = results{ii}.desired_obs(3);
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
% Get quantiles plus code uncertainty
post_cost_uq_cu = post_cost_uq + post_cost_sd;
post_cost_lq_cu = post_cost_lq - post_cost_sd;
post_defl_uq_cu = post_defl_uq + post_defl_sd;
post_defl_lq_cu = post_defl_lq - post_defl_sd;
post_rotn_uq_cu = post_rotn_uq + post_rotn_sd;
post_rotn_lq_cu = post_rotn_lq - post_rotn_sd;
% Get ylims for the two sets of plots
ylimrat=1.01;
ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];
ylim_rotn = [min(post_rotn_lq_cu)/ylimrat max(post_rotn_uq_cu)*ylimrat];

% Begin figures
h=figure('rend','painters','pos',[10 10 1200 400]);
x = 96:1:350; % x fills the cost domain
% Now begin plot 1/3
subplot(1,3,1)
% Get main curve
pdefl = pchip(cost,post_defl_mean,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(cost,post_defl_uq,x);
pdefllq = pchip(cost,post_defl_lq,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost,post_defl_mean,'or',...
    x,pdefl,'-r',...
    x,pdefluq,'-k',...
    x,pdefllq,'-k');
xl2=xlabel('Target cost');
ylabel('Deflection');
xlim([96,350]);
ylim(ylim_defl);

% Here's plot 2/3
subplot(1,3,2)
% Get main curve
protn = pchip(cost,post_rotn_mean,x);
% Get upper and lower 0.05 quantiles curves
protnuq = pchip(cost,post_rotn_uq,x);
protnlq = pchip(cost,post_rotn_lq,x);
f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost,post_rotn_mean,'or',...
    x,protn,'-r',...
    x,protnuq,'-k',...
    x,protnlq,'-k');
xl3=xlabel('Target cost');
ylabel('Rotation');
xlim([96,350]);
ylim(ylim_rotn);

% Here's plot 3/3
subplot(1,3,3)
% Get main curve
pcost = pchip(cost,post_cost_mean,x);
% Get upper and lower 0.05 quantiles curves
pcostuq = pchip(cost,post_cost_uq,x);
pcostlq = pchip(cost,post_cost_lq,x);
f=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost,post_cost_mean,'or',...
    x,pcost,'-r',...
    x,pcostuq,'-k',...
    x,pcostlq,'-k');
% Plot 2sd errbar
% errorbar(cost_lambda,post_cost_mean,post_cost_sd,'ob'); 
xl1=xlabel('Target cost');
ylabel('Observed cost');
xlim([96,350]);
ylim(ylim_cost);
%plot(x,x,'-k','LineWidth',2);

% Save the figure temporarily so we can mess with it later, because
% suptitle seems to mess things up somehow for making changes after calling
% it
savefig(h,'tempfig');

% Now add a main title and fix any infelicities
suptitle(['Performance metrics vs. (known) target cost,',...
    ' with ',num2str(cred_level),'% credible interval']); 
p = get(xl1,'position');
set(xl1,'position',p + [0 2.75 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.00125 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.0002 0])
figpos = get(h,'pos');

saveas(h,'FIG_costs.png');

% Now add in code uncertainty. That is, the above assumes that the GP
% emulator nails the FE code precisely. But of course the GP emulator has
% nonnegligible variance. That's the code uncertainty. So our confidence
% bands should reflect it. So we add it in here, by dropping the
% appropriate multiple of the sd from each lower quantile and adding it to
% each upper quantile.
% First, open the figure prior to calling suptitle.
h=openfig('tempfig');
subplot(1,3,1);
pdefluq_code_uncert = pchip(cost,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(cost,post_defl_lq_cu,x);
f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'b');
ff=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
xl2=xlabel('Target cost');
ylim(ylim_defl);

subplot(1,3,2);
protnuq_code_uncert = pchip(cost ,post_rotn_uq_cu,x);
protnlq_code_uncert = pchip(cost ,post_rotn_lq_cu,x);
f=fill([ x , fliplr(x) ], [protnuq_code_uncert,...
    fliplr(protnuq)],'b');
ff=fill([ x , fliplr(x) ], [protnlq_code_uncert,...
    fliplr(protnlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
xl3=xlabel('Target cost');
ylim(ylim_rotn);

subplot(1,3,3);
pcostuq_code_uncert = pchip(cost ,post_cost_uq_cu,x);
pcostlq_code_uncert = pchip(cost ,post_cost_lq_cu,x);
f=fill([ x , fliplr(x) ], [pcostuq_code_uncert,...
    fliplr(pcostuq)],'b');
ff=fill([ x , fliplr(x) ], [pcostlq_code_uncert,...
    fliplr(pcostlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
xl1=xlabel('Target cost');
ylim(ylim_cost);

% Now add a main title and fix any infelicities
suptitle(['Performance metrics vs. (known) target cost,',...
    ' with ',num2str(cred_level),'% credible interval ',...
    'including code uncertainty']); 
set(h,'pos',figpos); % Just so we can reuse the positioning code from above
p = get(xl1,'position');
set(xl1,'position',p + [0 2.75 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.00125 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.0002 0])

saveas(h,'FIG_costs_code_uncert.png');

