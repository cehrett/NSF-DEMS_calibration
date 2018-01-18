% Workspace: MCMC piecewise

clc; clear all; close all;

%% Get data
addpath('.\NSF DEMS\Phase 1\');
addpath('.\NSF DEMS\Phase 1\stored_data');

fprintf('Reading data from .xlsx...\n')
raw_dat = xlsread("fe_results.xlsx");
tdat = Tdat(raw_dat,3); % rescaling inputs, standardizing outputs
fprintf('done.\n\n')

%% Covariance parameter settings
% Found by optimization routine:
Rho_lam_optimum  = [0.655344235568109   0.931941001705886 ...
    0.960653924901867   0.991953924787049   0.017385994893994];
omega  = Rho_lam_optimum(1:2);
rho    = Rho_lam_optimum(3:4);
lambda = Rho_lam_optimum(5);

%% Set fake data values
defl_sd = .65/2;
defl_mean = .65;
rot_sd = .077/2;
rot_mean = .077;
cost_sd = 96/2;
cost_mean = 96;
% Set up inputs
temps = unique(raw_dat(:,1));
temps_std = (temps - tdat.input_mins(2))/tdat.input_ranges(2);
n =  length(temps);
ones_vec = ones(n,1);
cntrl_input = [0 * ones_vec temps_std ; .5 * ones_vec temps_std ; ...
    1 * ones_vec temps_std ];
% Set up outputs
% Note: since we're taking our fake data to be constant across temperature,
% there is a question of how to standardize it, since our original
% simulator data is not constant with respect to temperature, and was
% itself standardized at each temperature. Currently, to standardize the
% fake data, I use the mean of the by-temperature means, and the mean of
% the by-temperature sd's, and use these new means to standardize the fake
% data.
defl_mean_std = (defl_mean - mean(tdat.output_means(1,:)))/...
    mean(tdat.output_sds(1,:));
rot_mean_std = (rot_mean - mean(tdat.output_means(2,:)))/...
    mean(tdat.output_sds(2,:));
cost_mean_std = (cost_mean - mean(tdat.output_means(3,:)))/...
    mean(tdat.output_sds(3,:));
des_output = [ones_vec * defl_mean_std ; ones_vec * rot_mean_std ; ...
    ones_vec * cost_mean_std];

%% Settings for MCMC
M = 1e4 ; % Total number of draws
burn_in = ceil(M/5); % Total burn-in
init = rand(1,2);  % Initial calib parameter settings
%%%% Proposal density
prop_density = @(x,Sigma) randn * Sigma + x;
Sigma_v = 0.05;
Sigma_k = 0.05;
msg=0; % for console output
accepted_v = 0 ; % Use this for adaptive proposal distribution density
accepted_k = 0 ; % Use this for adaptive proposal distribution density
% The below vars will constrain uniform prior on VF, thickness
vf_min = tdat.input_mins(3); 
vf_range = tdat.input_ranges(3);
thk_min = tdat.input_mins(4);
thk_range = tdat.input_ranges(4);
init .* [vf_range thk_range] + [vf_min thk_min] % Take a look at the init
startplot = 1; % Start of samples to plot during MCMC

%% Create objects for data collection
samples = nan(M,2);
v = init(1);
k = init(2);
Sigma_xx  = gp_cov(omega,tdat.input(:,1:2),tdat.input(:,1:2),...
        rho,tdat.input(:,3:end),tdat.input(:,3:end),lambda);
Sigma_xx = Sigma_xx + eye(size(Sigma_xx))*1e-2;
inv_Sig_xx = inv(Sigma_xx);
% Get log likelihood of initial calibration value
pred_pts = [cntrl_input repmat([v k],size(cntrl_input,1),1)];
em=emulator(tdat.input,tdat.output,pred_pts,omega,rho,lambda,...
        Sigma_xx,inv_Sig_xx,false);
loglik_vk = logmvnpdf(des_output',em.mu',em.cov_gp);


% For debugging:
% sim_dat_input = tdat.input; sim_dat_output = tdat.output ; ii=1;
figure()
%% Begin MCMC
for ii = 1 : M
 
    % Propose new VF input
    vs = prop_density(v,Sigma_v);
    
    % Set up new prediction points based on xs
    pred_pts = [cntrl_input repmat([vs k],size(cntrl_input,1),1)];
    
    % Get (log) likelihood of vs
    if (vs>(.2-vf_min)/vf_range) && (vs< (.6-vf_min)/vf_range)
        em=emulator(tdat.input,tdat.output,pred_pts,omega,rho,lambda,...
            Sigma_xx,inv_Sig_xx,false);
        loglik_vs = logmvnpdf(des_output',em.mu',em.cov_gp) ;
    else loglik_vs=-Inf;
    end
    
    % Get acceptance ratio statistic
    log_alpha = loglik_vs - loglik_vk ;
    
    % Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted_v = accepted_v + 1;
    else 
        accept = 0;
    end
    
    if accept % Set up for next time
        loglik_vk = loglik_vs ;
        v = vs ;
    end
    
    % Propose new thickness input
    ks = prop_density(k,Sigma_k);
    
    % Set up new prediction points based on xs
    pred_pts = [cntrl_input repmat([v ks],size(cntrl_input,1),1)];
    
    % Get (log) likelihood of xs
    if (ks>(10-thk_min)/thk_range) &&(ks< (25-thk_min)/thk_range)
        em=emulator(tdat.input,tdat.output,pred_pts,omega,rho,lambda,...
            Sigma_xx,inv_Sig_xx,false);
        loglik_ks = logmvnpdf(des_output',em.mu',em.cov_gp) ;
    else loglik_ks = -Inf;
    end
    
    % Get acceptance ratio statistic
    log_alpha = loglik_ks - loglik_vk ;
    
    % Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted_k = accepted_k + 1;
    else 
        accept = 0;
    end
    
    if accept % Set up for next time
        loglik_vk = loglik_ks ;
        k = ks ;
    end
    
    % Recordkeeping
    samples(ii,:) = [v k];
    
    if mod(ii,10) == 0 % Output to console to let us know progress
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        subplot(1,2,1);
        plot(samples(startplot:end,1),'ko');
        subplot(1,2,2);
        plot(samples(startplot:end,2),'ko');
        drawnow
    end
    
    % This will effect adaptive proposal variance
    if mod(ii,100) == 0 && ii < burn_in
        fprintf(repmat('\b',1,msg));
        if accepted_v < 40 
            Sigma_v = Sigma_v * 0.75;
            fprintf('VF proposal Variance reduced to %g\n',Sigma_v);
        end
        if accepted_v > 50
            Sigma_v = Sigma_v * 1.25;
            fprintf('VF proposal Variance increased to %g\n',Sigma_v);
        end
        if accepted_k < 40 
            Sigma_k = Sigma_k * 0.75;
            fprintf('Thickness proposal Variance reduced to %g\n',Sigma_k);
        end
        if accepted_k > 50
            Sigma_k = Sigma_k * 1.25;
            fprintf('Thickness proposal Variance increased to %g\n',Sigma_k);
        end
        msg = fprintf('Completed: %g/%g\n',ii,M);
        accepted_v = 0;
        accepted_k = 0;
    end
    
    % Stop plotting burn_in
    if ii > burn_in
        startplot=burn_in;
    end
    
end
% 
% samples
% mean(samples(burn_in:end,:))
% var(samples(burn_in:end,:))
% figure()
% plot(samples(burn_in:end,1),'bo')
% figure()
% plot(samples(burn_in:end,2))

pbi_samples = samples(burn_in:end,:);
thinning = 1:5:length(pbi_samples);
pbi_samples(thinning,:)
post_mean_thinned = mean(pbi_samples(thinning,:))
post_var_thinned  = var(pbi_samples(thinning,:))
figure()
plot(pbi_samples(thinning,1));
figure()
plot(pbi_samples(thinning,2));

Sigma = [Sigma_v Sigma_k];
save('.\NSF DEMS\Phase 1\samples6','samples');
save('.\NSF DEMS\Phase 1\init6','init');
save('.\NSF DEMS\Phase 1\Sigma6','Sigma');
Sigma

figure()
o_samples = pbi_samples;
o_samples(:,1) = o_samples(:,1)*vf_range + vf_min;
o_samples(:,2) = o_samples(:,2)*tdat.input_ranges(4) + tdat.input_mins(4);
subplot(1,2,1);
plot(o_samples(:,1),'ko');
title('Volume fraction');
subplot(1,2,2);
plot(o_samples(:,2),'ko');
title('Thickness');

% Load all three saved samples and plot them
load('.\NSF DEMS\Phase 1\samples1');
samples(:,1) = samples(:,1) * vf_range + vf_min;
samples(:,2) = samples(:,2) * thk_range + thk_min;
samples1=samples;
load('.\NSF DEMS\Phase 1\samples2');
samples(:,1) = samples(:,1) * vf_range + vf_min;
samples(:,2) = samples(:,2) * thk_range + thk_min;
samples2=samples;
load('.\NSF DEMS\Phase 1\samples3');
samples(:,1) = samples(:,1) * vf_range + vf_min;
samples(:,2) = samples(:,2) * thk_range + thk_min;
samples3=samples;
figure()
subplot(1,2,1);
plot(samples1(burn_in:end,1),'ko');
subplot(1,2,2);
plot(samples1(burn_in:end,2),'ko');
figure()
subplot(1,2,1);
plot(samples(burn_in:end,1),'ko');
subplot(1,2,2);
plot(samples(burn_in:end,2),'ko');
figure()
subplot(1,2,1);
plot(samples(burn_in:end,1),'ko');
subplot(1,2,2);
plot(samples(burn_in:end,2),'ko');

