% Workspace: MCMC
% This version has a joint proposal density, and uses Z=[Y;eta] in the
% MCMC routine, where Y is transformed version of W; W Exp(eta), Y normal.

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
x = [0 * ones_vec temps_std ; .5 * ones_vec temps_std ; ...
    1 * ones_vec temps_std ];
% Set up outputs
% Note: since we're taking our fake data to be constant across temperature,
% there is a question of how to standardize it, since our original
% simulator data is not constant with respect to temperature, and was
% itself standardized at each temperature. Currently, to standardize the
% fake data, I use the mean of the by-temperature means, and the mean of
% the by-temperature sd's, and use these new means to standardize the fake
% data. Ultimately this should not matter much, because the means and sd's
% do not vary much across temperature.
defl_mean_std = (defl_mean - mean(tdat.output_means(1,:)))/...
    mean(tdat.output_sds(1,:));
defl_0_std = (- mean(tdat.output_means(1,:)))/...
    mean(tdat.output_sds(1,:));
defl_sd_std = (defl_mean_std - defl_0_std) / 2 ; 
rot_mean_std = (rot_mean - mean(tdat.output_means(2,:)))/...
    mean(tdat.output_sds(2,:));
rot_0_std = (- mean(tdat.output_means(2,:)))/...
    mean(tdat.output_sds(2,:));
rot_sd_std = (rot_mean_std - rot_0_std) / 2 ; 
cost_mean_std = (cost_mean - mean(tdat.output_means(3,:)))/...
    mean(tdat.output_sds(3,:));
cost_0_std = (- mean(tdat.output_means(3,:)))/...
    mean(tdat.output_sds(3,:));
cost_sd_std = (cost_mean_std - cost_0_std) / 2 ; 
des_output = [ones_vec * defl_mean_std ; ones_vec * rot_mean_std ; ...
    ones_vec * cost_mean_std];
W = 1e-4*ones(size(des_output)); % We'll transform this to normal; for now, we consider W_i 
                % to be distrib. as Exp(eta(x_i,theta)).

%% Settings for MCMC
M = 1e3 ; % Total number of draws
xts = tdat.input ; % Sim input
n = size(W,1); % Total number of "observations"
m = size(xts,1); % Total number of simulation runs
burn_in = floor(M/5); % Total burn-in
init = rand(1,2)  % Initial calib parameter settings
%%%% OPTIONS for proposal density: Beta (bdd), Normal on logit tran (bdd)...
%%%% or just Normal (unbdd)
prop_dens_ind = 1;
prop_density = @(x,n) (betarnd(n*x,(n*x-n*x.*x)./x) ); % Beta proposal dens.
% Logit transformations for logit trans normal proposal
prop_dens_ind = 2;
logit_trans = @(x) log(x ./ (1-x));
logit_rev_trans = @(x) exp(x) ./ (1+exp(x));
prop_density = @(x,Sigma) logit_rev_trans(mvnrnd(logit_trans(x),Sigma));
Sigma = [.1 0 ; 0 .1];
prop_dens_ind = 3;
prop_density = @(x,Sigma) (mvnrnd(x,Sigma)); % Normal
Sigma = [.01 0 ; 0 .01];
msg=0; % for console output
accepted = 0 ; % Use this for adaptive proposal distribution density
% The below vars will constrain uniform prior on VF, thickness
vf_min = tdat.input_mins(3); 
vf_range = tdat.input_ranges(3);
thk_min = tdat.input_mins(4);
thk_range = tdat.input_ranges(4);
startplot = 1; % Start of samples to plot during MCMC

%% Create objects for data collection
samples = nan(M,2);
v = init(1);
k = init(2);
theta = init;
% Transform W into Y~N(eta,sigma^2)
% First we need (stand-ins for) eta(x_i,theta) values.
% em = emulator(xts,tdat.output,[x repmat(theta,size(x,1),1)],omega,rho,...
%    lambda,0,0,true);
% eta = em.mu;
% return eta to original scale
means=[repmat(mean(tdat.output_means(1,:)),n/3,1); ...
    repmat(mean(tdat.output_means(2,:)),n/3,1);    ...
    repmat(mean(tdat.output_means(3,:)),n/3,1)];
% means(1:n/3) = mean(tdat.output_means(1,:));
% means(1+n/3:2*n/3) = mean(tdat.output_means(2,:));
% means(1+2*n/3:end) = mean(tdat.output_means(3,:));
sds=[repmat(mean(tdat.output_sds(1,:)),n/3,1); ...
    repmat(mean(tdat.output_sds(2,:)),n/3,1);    ...
    repmat(mean(tdat.output_sds(3,:)),n/3,1)];
% eta_or = eta .* sds + means;
eta_or = [.65 ;.077 ;96 ] ; eta_or=repelem(eta_or,21) ;
% Next we use these to transform W
% First, we get Wq = F(W), where F is the cdf of W.
Wq = 1 - exp(-W./eta_or) ;
% Next, we get the inverse normal(eta_or,Sigma) quantiles of Wq. To do
% that, we need a vector of the diagonal of Sigma (which is a diag matrix)
Sigma_diag = [repmat(defl_sd_std,n/3,1); repmat(rot_sd_std,n/3,1); ...
    repmat(cost_sd_std,n/3,1)];
Y = norminv(Wq,eta_or,Sigma_diag);
Y = (Y-means)./Sigma_diag ;
z = [Y ; tdat.output ] ; % Vector of observations and sim output
inputs = [ x repmat(theta,size(x,1),1) ; xts ] ;
% Get Sigma_eta
Sigma_eta_yy = gp_cov(omega,inputs(1:n,1:2),inputs(1:n,1:2),rho,...
    inputs(1:n,3:4),inputs(1:n,3:4),lambda,true);
Sigma_eta_xy = gp_cov(omega,inputs(n+1:n+m,1:2),inputs(1:n,1:2),rho,...
    inputs(n+1:n+m,3:4),inputs(1:n,3:4),lambda,true);
Sigma_eta_yx = Sigma_eta_xy';
Sigma_eta_xx = gp_cov(omega,inputs(n+1:n+m,1:2),inputs(n+1:n+m,1:2),...
    rho,inputs(n+1:n+m,3:4),inputs(n+1:n+m,3:4),lambda,true);
Sigma_eta = [Sigma_eta_yy Sigma_eta_yx ; Sigma_eta_xy Sigma_eta_xx ] ;
% Now build Sigma_y
var_y = [defl_sd_std rot_sd_std cost_sd_std].^2;
var_y = repelem(var_y,n/3);
Sigma_y = diag(var_y);
% Now we can get Sigma_z
Sigma_z = Sigma_eta + padarray(Sigma_y,[m m],'post') ;
% The following 2 lines are for numerical stability
nugsize = 1e-4 ;
Sigma_z = Sigma_z + eye(size(Sigma_z))*nugsize;
% Get the initial loglikelihood.
loglik_theta = logmvnpdf(z',0,Sigma_z);
vf_outofbounds = 0; % These will help keep track for tuning adaptive var
thk_outofbounds =0;
reject = false;

% For debugging:
% sim_dat_input = tdat.input; sim_dat_output = tdat.output ; ii=1;
figure()
%% Begin MCMC
for ii = 1 : M
 
    % Propose new calib input
    theta_s = prop_density(theta,Sigma);
    
    % Get (log) likelihood of new input
    if theta_s(1) < (.2-vf_min)/vf_range || theta_s(1) > ...
            (.6-vf_min)/vf_range
        vf_outofbounds = vf_outofbounds + 1 ;
        reject = true ;
    end
    if theta_s(2) < (10-thk_min)/thk_range || theta_s(2) > ...
            (25-thk_min)/thk_range
        thk_outofbounds = thk_outofbounds + 1;
        reject = true;
    end
    if reject == false
        inputs(1:n,3:4) = repmat(theta_s,n,1) ;
        % Get new Sigma_eta_yy 
        Sigma_eta_yy = gp_cov(omega,inputs(1:n,1:2),inputs(1:n,1:2),rho,...
            inputs(1:n,3:4),inputs(1:n,3:4),lambda,false);
        % Get new Sigma_eta_xy, and hence Sigma_eta_yx
        Sigma_eta_xy = gp_cov(omega,inputs(n+1:n+m,1:2),inputs(1:n,1:2),...
            rho,inputs(n+1:n+m,3:4),inputs(1:n,3:4),lambda,true);
        Sigma_eta_yx = Sigma_eta_xy';
        % Combine these to get new Sigma_z
        Sigma_z = [ Sigma_eta_yy+Sigma_y  Sigma_eta_yx ; Sigma_eta_xy ...
            Sigma_eta_xx ] ;
        %Sigma_eta = gp_cov(omega,inputs(:,1:2),inputs(:,1:2),rho,inputs(:,3:4),...
        %inputs(:,3:4),lambda);
        %Sigma_z = Sigma_eta + padarray(Sigma_y,[m m],'post') ;
        % The following line is for numerical stability
        Sigma_z = Sigma_z + eye(size(Sigma_z))*nugsize;
        loglik_theta_s = logmvnpdf(z',0,Sigma_z) ;
    else
        loglik_theta_s = -Inf;
        reject = false;
    end
    
    % Get ratio of (log) likelihoods of theta|theta_s, theta_s|theta under
    % proposal density
    % Note that the below assumes the prop density is (independent) 
    % normal distributions of logit-transformed calibration input variables
    if prop_dens_ind == 2
        loglik_prop_ratio = sum(log(theta_s .* (1-theta_s)) - ...
            log(theta .* (1-theta)));
    else loglik_prop_ratio = 0;
    end
    
    % Get acceptance ratio statistic
    log_alpha = loglik_theta_s - loglik_theta + loglik_prop_ratio;
    
    % Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted = accepted + 1;
    else 
        accept = 0;
    end
    
    if accept % Set up for next time
        loglik_theta = loglik_theta_s ;
        theta = theta_s ;
    end
    
    % Recordkeeping
    samples(ii,:) = theta;
    
    if mod(ii,10) == 0 % Output to console to let us know progress
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        subplot(1,2,1);
        plot(samples(startplot:end,1),'ko');
        subplot(1,2,2);
        plot(samples(startplot:end,2),'ko');
        drawnow
    end
    
    if mod(ii,100) == 0 && ii < burn_in % Tune adaptive proposal variance
        if vf_outofbounds >= 40
            Sigma(1,1) = .75 * Sigma(1,1);
            fprintf(repmat('\b',1,msg));
            fprintf('VF proposal variance reduced to %g\n',Sigma(1,1));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if thk_outofbounds >= 40
            Sigma(2,2) = .75 * Sigma(2,2);
            fprintf(repmat('\b',1,msg));
            fprintf('Thk proposal variance reduced to %g\n',Sigma(2,2));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if vf_outofbounds < 40 && thk_outofbounds < 40 && accepted < 30 
            Sigma = Sigma * 0.75;
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances reduced to %g,%g\n',diag(Sigma));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted > 50
            Sigma = Sigma * 1.25;
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances increased to %g\n',Sigma(1,1));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        accepted        = 0;
        vf_outofbounds  = 0;
        thk_outofbounds = 0;
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

save('.\NSF DEMS\Phase 1\expsamples1','samples');
save('.\NSF DEMS\Phase 1\expinit1','init');
save('.\NSF DEMS\Phase 1\expSigma1','Sigma');
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

