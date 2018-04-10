% Workspace for grid investigation of Cost_lambda

clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

%% User defined values
m = 10; % Grid size
M = 6e3; % MCMC length
desired_obs = [ 0 0 ];
which_outputs = [ 1 1 0 ] ; %Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
settings.burn_in=2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop through grid of Cost_lambda vals
%gridvals = exp(linspace(0,log(100),m));
gridvals = linspace(0,100,m); % Get grid values
results = cell(m,1); % This will store results
for ii = 1:4%1:m
    
    settings.Cost_lambda = gridvals(ii);
    
    % Joint prop for theta, joint prop for obs var, prior on obs var
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    
    % Collect results
    % trivar_output_settings is useful for getting output at post mean
    trivar_output_settings = MCMC_settings(M,[0 0 0],[1 1 1]);
    post_mean_out = em_out(samples,trivar_output_settings);
    result = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
        'init',samples(1,:),'desired_obs',desired_obs,...
        'sigma2_prior',settings.log_sigma2_prior,...
        'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
        'proposal',settings.proposal,'nugsize',settings.nugsize,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,'which_outputs',which_outputs,...
        'Cost_lambda',settings.Cost_lambda,...
        'log_theta_prior',settings.log_theta_prior,...
        'log_sig_mh_correction',settings.log_sig_mh_correction,...
        'log_mh_correction',settings.log_mh_correction,...
        'burn_in',settings.burn_in);
    
    results{ii+1} = result;
    
    save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
     'results_Cost_lambda_grid_exploration'],...
     'results');
    
%     save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
%     'results_Cost_lambda_grid_exploration'],...
%     'results');
    
end

% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
% 'results_Cost_lambda_grid_exploration'],...
% 'results');

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
'results_Cost_lambda_grid_exploration'],...
'results');

%% Figures
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
 'results_Cost_lambda_grid_exploration'],...
 'results');
% Collect Cost_lambdas, posterior costs, defl, rot
m=size(results,1);
cost_lambda = zeros(m,1);
pmo = zeros(m,3);
for ii = 1:m
    pmo(ii,:) = results{ii}.post_mean_out;
    cost_lambda(ii) = results{ii}.Cost_lambda;
end
post_cost = pmo(:,3);
post_defl = pmo(:,1);
post_rotn = pmo(:,2);
% Begin figures
h=figure('rend','painters','pos',[10 10 1200 400]);
x = 0:1:100;
subplot(1,3,1)
pcost = pchip(cost_lambda,post_cost,x);
plot(cost_lambda,post_cost,'o',x,pcost,'-');
xl1=xlabel('\lambda_c_o_s_t');
ylabel('Cost');

subplot(1,3,2)
pdefl = pchip(cost_lambda,post_defl,x);
plot(cost_lambda,post_defl,'o',x,pdefl,'-');
xl2=xlabel('\lambda_c_o_s_t');
ylabel('Deflection');

subplot(1,3,3)
protn = pchip(cost_lambda,post_rotn,x);
plot(cost_lambda,post_rotn,'o',x,protn,'-');
xl3=xlabel('\lambda_c_o_s_t');
ylabel('Rotation');


suptitle('Performance metrics vs. \lambda_c_o_s_t'); 
p = get(xl1,'position');
set(xl1,'position',p + [0 6 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.003 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.0005 0])

saveas(h,'FIG_cost_lambda.png');

%% Get confidence intervals at each Cost_lambda
m=size(results,1);

intervals = 
