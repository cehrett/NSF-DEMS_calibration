function settings = MCMC_dual_calib_settings(...
    sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,...
    des_x,des_y,...
    varargin)
% settings = MCMC_settings(sim_x,sim_t1,sim_t2,sim_y,...
%                          obs_x,obs_t2,obs_y,des_x,des_y,varargin)
% 
% Parameter/Value pairs:
% 'M':
%   Number of samples to draw in the MCMC, includes burn-in. Default: 1e4.
% 'burn_in':
%   Proportion of samples to treat as burn-in. Default: 1/5.
% 'min_x':
%   Vector giving the minimum value of each element of the control input.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_x and the observation input obs_x.
% 'range_x':
%   Vector giving the ranges of each element of the control input.
%   By default, this is taken to be the range of the simulator input
%   provided as sim_x and the observation input obs_x.
% 'min_t1':
%   Vector giving the minimum value of each element of theta1 parameter.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_t1.
% 'range_t1':
%   Vector giving the range of each element of theta1 parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t1.
% 'dim_t1':
%   Dimensionality of theta1 calibration parameter. By default, this is
%   inferred from sim_t1; or, if an emulator is not used and sim_t1 is not
%   supplied, it is assumed to be 1.
% 'min_t2':
%   Vector giving the minimum value of each element of theta2 parameter.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'range_t2':
%   Vector giving the range of each element of theta2 parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'ObsVar':
%   Scalar value which gives the variance of the iid true observation
%   error. Default: 0.05.
% 'EmulatorCovHypers':
%   Vector of hyperparameters for the power product exponential covariance
%   function of the GP emulator. Default:
%               [ 0.5 0.5 0.5 Inf ]. 
%   If set to 0, then ML estimation via gradient descent is used to find 
%   appropriate values.
% 'DiscrepancyCovHypers':
%   Vector of hyperparameters for the power product exponential covariance
%   function of the GP discrepancy function (for the discrepancy between
%   the simulator and the true observations). The vector should consist of
%   the rho values (rho_i = exp(-1/(2*L_i^2) where L_i is the lengthscale
%   of the covariance in the ith dimension of inputs) followed by lambda,
%   the marginal precision. Default:
%               [ 0.5 0.5 Inf ]. 
%   If set to 0, then ML estimation via gradient descent is used to find 
%   appropriate values.
% 'doplot':
%   true       - (Default) Update scatterplot every ten draws.
%   false      - No plots during MCMC.


%% Parse inputs
p = inputParser;
p.addRequired('sim_x',@ismatrix);
p.addRequired('sim_t1',@ismatrix);
p.addRequired('sim_t2',@ismatrix);
p.addRequired('sim_y',@ismatrix);
p.addRequired('obs_x',@ismatrix);
p.addRequired('obs_t2',@ismatrix);
p.addRequired('obs_y',@ismatrix);
p.addRequired('des_x',@ismatrix);
p.addRequired('des_y',@ismatrix);
p.addParameter('M',1e4,@isscalar);
p.addParameter('burn_in',1/5,@(x) x>0 && x<1);
p.addParameter('min_x','Default',@isnumeric);
p.addParameter('range_x','Default',@isnumeric);
p.addParameter('min_t1','Default',@isnumeric);
p.addParameter('range_t1','Default',@isnumeric);
p.addParameter('dim_t1',size(sim_t1,2),@isscalar);
p.addParameter('min_t2','Default',@isnumeric);
p.addParameter('range_t2','Default',@isnumeric);
p.addParameter('ObsVar',0.05,@isscalar);
p.addParameter('EmulatorCovHypers',[0.5 0.5 0.5 Inf],@ismatrix);
p.addParameter('DiscrepancyCovHypers',[0.5 0.5 Inf],@ismatrix);
p.addParameter('doplot',true,@islogical);
p.parse(sim_x,sim_t1,sim_t2,sim_y,obs_x,obs_t2,obs_y,des_x,des_y,...
    varargin{:});


%% Collect inputs
M = p.Results.M;
burn_in = floor(p.Results.burn_in*M);
min_x = p.Results.min_x;
range_x = p.Results.range_x;
min_t1 = p.Results.min_t1;
range_t1 = p.Results.range_t1;
dim_t1 = p.Results.dim_t1; if isequal(sim_t1,[]), dim_t1=1; end
min_t2 = p.Results.min_t2;
range_t2 = p.Results.range_t2;
ObsVar = p.Results.ObsVar;
EmulatorCovHypers = p.Results.EmulatorCovHypers;
DiscrepancyCovHypers = p.Results.DiscrepancyCovHypers;
doplot = p.Results.doplot;


%% Some useful definitions
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x) ./ (1+exp(x));


%% Normalize inputs
x = [sim_x ; obs_x] ;
if min_x == 'Default', min_x = min(x) ; end
if range_x == 'Default', range_x = range(x) ; end
sim_x_01 = (sim_x - min_x) ./ range_x ; 
obs_x_01 = (obs_x - min_x) ./ range_x ;
des_x_01 = (des_x - min_x) ./ range_x ;

if min_t1 == 'Default', min_t1 = min(sim_t1) ; end
if range_t1 == 'Default', range_t1 = range(sim_t1) ; end
sim_t1_01 = (sim_t1 - min_t1) ./ range_t1 ;

t2 = [sim_t2 ; obs_t2 ] ;
if min_t2 == 'Default', min_t2 = min(t2) ; end
if range_t2 == 'Default', range_t2 = range(t2) ; end
sim_t2_01 = (sim_t2 - min_t2) ./ range_t2 ; 
obs_t2_01 = (obs_t2 - min_t2) ./ range_t2 ;


%% Standardize outputs
y = [sim_y ; obs_y] ; 
mean_y = mean(y) ; 
std_y = std(y) ; 
sim_y_std = (sim_y - mean_y) ./ std_y ; 
obs_y_std = (obs_y - mean_y) ./ std_y ;
des_y_std = (des_y - mean_y) ./ std_y ; 


%% Set emulator mean
% This is currently set so that an emulator is, strictly speaking, not
% used. That is, the emulator mean function is set to just be the objective
% function itself, with the expectation that this will be paired with a 0
% covariance function.
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,min_x,range_x,b,min_t1,range_t1,c,min_t2,range_t2,...
    mean_y,std_y);
%mean_sim = @(a,b,c) zeros(size(a,1)); % Emulator mean


%% Make true observation error covariance matrix
num_obs = size(obs_y(:),1) ; 
obs_cov_mat = ObsVar * eye(num_obs) ; 


%% Set theta priors
log_theta1_prior = @(t) 0 ; % Uniform prior
log_theta2_prior = @(t) 0 ; % Uniform prior


%% Set theta proposal distributions and initial values
% We'll draw logit-transformed theta from a normal centered at the
% logit-transformed previous draw.
theta1_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
theta2_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
% Since this proposal is not symmetric, we need to use full
% Metropolis-Hastings rather than just Metropolis. So here is the log MH
% correction for the lack of symmetry.
theta1_prop_log_mh_correction = ...
    @(t_s,t) sum(log(t_s)+log(1-t_s)-log(t)-log(1-t));
theta2_prop_log_mh_correction = ...
    @(t_s,t) sum(log(t_s)+log(1-t_s)-log(t)-log(1-t));
% Set initial values and initial covariance matrices for proposals
theta1_init = rand(dim_t1,1);
theta2_init = rand(dim_t1,1);
Sigma_theta1 = eye(size(theta1_init,1));
Sigma_theta2 = eye(size(theta2_init,1));


%% Set rho and lambda prior distributions
log_des_discrep_rho_prior  = @(r) sum(log( betapdf(r,1,0.6) ));
log_des_discrep_lambda_prior = @(ld) log( gampdf(ld,5,5) );


%% Set rho and lambda proposal distributions
% We'll draw logit-transformed rho and lambda from normals centered at the
% log-transformed previous draw.
rho_proposal = @(r,S) logit_inv(mvnrnd(logit(r),S)); 
lambda_proposal = @(lam,S) exp(mvnrnd(log(lam),S)); 
% Since this proposal is not symmetric, we need to use full
% Metropolis-Hastings rather than just Metropolis. So here is the log MH
% correction for the lack of symmetry.
rho_prop_log_mh_correction = ...
    @(r_s,r) sum(log(r_s)+log(1-r_s)-log(r)-log(1-r));
lambda_prop_log_mh_correction = ...
    @(lam_s,lam) sum(log(lam_s)-log(lam));
% Set initial values and initial covariance matrices for proposals
rho_init = rand(size(des_x,2),1);
lambda_init = gamrnd(1,1);
Sigma_rho = eye(size(rho_init,1));
Sigma_lambda = 1;


%% Pack up the settings structure
settings = struct(...
    'M',M,...
    'burn_in',burn_in,...
    'sim_x',sim_x_01,...
    'sim_t1',sim_t1_01,...
    'sim_t2',sim_t2_01,...
    'sim_y',sim_y_std,...
    'obs_x',obs_x_01,...
    'obs_t2',obs_t2_01,...
    'obs_y',obs_y_std,...
    'des_x',des_x_01,...
    'des_y',des_y_std,...
    'min_x',min_x,...
    'range_x',range_x,...
    'min_t1',min_t1,...
    'range_t1',range_t1,...
    'min_t2',min_t2,...
    'range_t2',range_t2,...
    'mean_y',mean_y,...
    'std_y',std_y,...
    'obs_cov_mat',obs_cov_mat,...
    'mean_sim',mean_sim,...
    'emulator_rho',EmulatorCovHypers(1:end-1),...
    'emulator_lambda',EmulatorCovHypers(end),...
    'obs_discrep_rho',DiscrepancyCovHypers(1:end-1),...
    'obs_discrep_lambda',DiscrepancyCovHypers(end),...
    'rho_proposal',rho_proposal,...
    'lambda_proposal',lambda_proposal,...
    'rho_prop_log_mh_correction',rho_prop_log_mh_correction,...
    'lambda_prop_log_mh_correction',lambda_prop_log_mh_correction,...
    'des_discrep_rho_init',rho_init,...
    'des_discrep_lambda_init',lambda_init,...
    'rho_prop_cov',Sigma_rho,...
    'lambda_prop_cov',Sigma_lambda,...
    'des_discrep_log_rho_prior',log_des_discrep_rho_prior,...
    'des_discrep_log_lambda_prior',log_des_discrep_lambda_prior,...
    'theta1_init',theta1_init,...
    'theta2_init',theta2_init,...
    'theta1_proposal',theta1_proposal,...
    'theta2_proposal',theta2_proposal,...
    'theta1_prop_log_mh_correction',theta1_prop_log_mh_correction,...
    'theta2_prop_log_mh_correction',theta2_prop_log_mh_correction,...
    'theta1_prop_cov',Sigma_theta1,...
    'theta2_prop_cov',Sigma_theta2,...
    'log_theta1_prior',log_theta1_prior,...
    'log_theta2_prior',log_theta2_prior,...
    'doplot',doplot);

end