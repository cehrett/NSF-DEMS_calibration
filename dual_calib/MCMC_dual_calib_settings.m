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
%   Vector giving the minimum value of each element of the parameter to be
%   calibrated. By default, this is taken to be the minimum of the 
%   simulator input provided as sim_t1.
% 'range_t1':
%   Vector giving the range of each element of calibration parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t1.
% 'dim_x':
%   Dimensionality of x control input. By default, this is
%   inferred from the x inputs.
% 'dim_t1':
%   Dimensionality of theta1 calibration parameter. By default, this is
%   inferred from the t1 inputs.
% 'dim_t2':
%   Dimensionality of theta2 design parameter. By default, this is
%   inferred from the t2 inputs.
% 'min_t2':
%   Vector giving the minimum value of each element of design parameter.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'range_t2':
%   Vector giving the range of each element of design parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'mean_y':
%   Scalar giving the mean value of the output y. By default, this is taken
%   to be the mean of the simulator and observational data sim_y and obs_y.
% 'std_y':
%   Scalar giving the standard deviation of the output y. By default, this
%   is taken to be the standard deviation of the simulator and
%   observational data sim_y and obs_y.
% 'obs_var':
%   Gives the covariance matrix of the observation error. Inputting a
%   scalar value here results in a covariance matrix that is the identity
%   multiplied by that scalar value. Default: scalar 0.05.
% 'des_var':
%   Gives the covariance matrix of the target outcomes' error. Inputting a
%   scalar value here results in a covariance matrix that is the identity
%   multiplied by that scalar value. Default: scalar 0.05.
% 'additional_discrep_mean':
%   Optional additional discrepancy mean. This is primarily intended to
%   be used in the case when CTO is completed after KOH; this additional
%   discrepancy term allows us to include in CTO information about the
%   observation discrepancy learned during KOH. Default 0.
% 'additional_discrep_cov':
%   Optional additional discrepancy covariance. This is primarily intended 
%   to be used in the case when CTO is completed after KOH; this additional
%   discrepancy term allows us to include in CTO information about the
%   observation discrepancy learned during KOH. Default 0.
% 'emulator':
%   Determines whether an emulator is used. If not, the emulator mean 
%   function is set to just be the objective function itself, and pairs 
%   this with a 0 covariance fn (by setting marginal precision to Inf).
% 'EmulatorMean':
%   Mean function of the emulator. Default: constant 0 (if an emulator is
%   used), example objective function (if emulator==false).
% 'EmulatorCovHypers':
%   Vector of hyperparameters for the power product exponential covariance
%   function of the GP emulator. Default:
%               [0.992943679103582   0.785517245465518   ...
%                0.077856518100309   0.083606519464691 ]. 
% 'obs_discrep':
%   Boolean value which tells whether or not to include a discrepancy term
%   for the real observations. Default: true.
% 'des_discrep':
%   Boolean value which tells whether or not to include a discrepancy term
%   for the desired observations. Default: true.
% 'obs_discrep_mean':
%   Function of x and theta2 which gives the prior mean on the observation
%   discrepancy function. Default: constant 0.
% 'obs_rho_lambda':
%   Allows user to specify a value for obs_rho,obs_lambda By default, these
%   are estimated via MCMC.
% 'modular':
%   Boolean value which tells whether or not to use modularized version of
%   the model, to protect the traditional calibration from being influenced
%   by the desired observations. Default: false.
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
p.addParameter('dim_x','Default',@isscalar);
p.addParameter('dim_t1','Default',@isscalar);
p.addParameter('dim_t2','Default',@isscalar);
p.addParameter('min_t2','Default',@isnumeric);
p.addParameter('range_t2','Default',@isnumeric);
p.addParameter('mean_y','Default',@isscalar);
p.addParameter('std_y','Default',@isscalar);
p.addParameter('obs_var',0.05,@isscalar);
p.addParameter('des_var',0,@isscalar);
p.addParameter('additional_discrep_mean',@(x)0,...
    @(h)isa(h,'function_handle'));
p.addParameter('additional_discrep_cov',@(x)0,...
    @(h)isa(h,'function_handle'));
p.addParameter('emulator',true,@islogical);
p.addParameter('EmulatorMean','Default',@(h)isa(h,'function_handle'));
p.addParameter('EmulatorCovHypers','Default',@ismatrix);
p.addParameter('obs_discrep',true,@islogical);
p.addParameter('des_discrep',true,@islogical);
p.addParameter('obs_discrep_mean','Default',@(h)isa(h,'function_handle'));
p.addParameter('obs_rho_lambda','Default',@isnumeric);
p.addParameter('modular',false,@islogical);
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
dim_x = p.Results.dim_x; 
dim_t1 = p.Results.dim_t1; 
dim_t2 = p.Results.dim_t2;
min_t2 = p.Results.min_t2;
range_t2 = p.Results.range_t2;
mean_y = p.Results.mean_y;
std_y = p.Results.std_y;
obs_var = p.Results.obs_var;
des_var = p.Results.des_var;
additional_discrep_cov = p.Results.additional_discrep_cov;
additional_discrep_mean = p.Results.additional_discrep_mean;
emulator = p.Results.emulator;
EmulatorMean = p.Results.EmulatorMean;
EmulatorCovHypers = p.Results.EmulatorCovHypers;
obs_discrep = p.Results.obs_discrep;
des_discrep = p.Results.des_discrep;
obs_discrep_mean = p.Results.obs_discrep_mean;
obs_rho_lambda = p.Results.obs_rho_lambda;
modular = p.Results.modular;
doplot = p.Results.doplot;

%% Todo
if ~(isequal(obs_rho_lambda,'Default'))
    error('Error: obs_rho_lambda input parameter not implemented yet!');
end


%% Some useful definitions
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x) ./ (1+exp(x));


%% Normalize inputs
x = [sim_x ; obs_x] ;
if isequal(min_x,'Default'), min_x = min(x) ; end
if isequal(range_x,'Default'), range_x = range(x) ; end
if sum(size(sim_x))>0, sim_x_01 = (sim_x - min_x) ./ range_x ;
else, sim_x_01 = sim_x ; end
obs_x_01 = (obs_x - min_x) ./ range_x ;
if sum(size(des_x))>0, des_x_01 = (des_x - min_x) ./ range_x ;
else, des_x_01 = des_x ; end

if min_t1 == 'Default', min_t1 = min(sim_t1) ; end
if range_t1 == 'Default', range_t1 = range(sim_t1) ; end
sim_t1_01 = (sim_t1 - min_t1) ./ range_t1 ;

t2 = [sim_t2 ; obs_t2 ] ;
if isequal(min_t2,'Default'), min_t2 = min(t2) ; end
if isequal(range_t2,'Default'), range_t2 = range(t2) ; end
sim_t2_01 = (sim_t2 - min_t2) ./ range_t2 ; 
obs_t2_01 = (obs_t2 - min_t2) ./ range_t2 ;


%% Infer some useful values
if isequal(dim_x,'Default')
    dim_x = numel(min_x);
end
if isequal(dim_t1,'Default')
    dim_t1 = numel(min_t1);
end
if isequal(dim_t2,'Default')
    dim_t2 = numel(min_t2);
end



%% Standardize outputs
if mean_y == 'Default', mean_y = mean(sim_y) ; end
if std_y == 'Default', std_y = std(sim_y) ; end
sim_y_std = (sim_y - mean_y) ./ std_y ; 
obs_y_std = (obs_y - mean_y) ./ std_y ;
des_y_std = (des_y - mean_y) ./ std_y ; 


%% Set emulator mean and covariance hyperparameters
% Determines whether an emulator is used. If not, the emulator 
% mean function is set to just be the objective
% function itself, and pairs this with a 0
% covariance function (by setting marginal precision to Inf).
if emulator
    if isequal(EmulatorMean,'Default')
        mean_sim = @(a,b,c) zeros(size(a,1),1); % Emulator mean
    else
        mean_sim = EmulatorMean;
    end
    if isequal(EmulatorCovHypers,'Default')
        % Define function for minimization
        f = @(rl) ...
            -logmvnpdf(((sim_y-mean_y)./std_y)',...
            mean_sim(sim_x_01,sim_t1_01,sim_t2_01)',...
            gp_cov(rl(1:(end-1)),...
            [sim_x_01 sim_t1_01 sim_t2_01],...
            [sim_x_01 sim_t1_01 sim_t2_01],...
            rl(end),false) + ...
            1e-4*eye(size(sim_x_01,1)));
        % Perform minimization
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [zeros(1,size([sim_x sim_t1 sim_t2],2)) 0];
        ub = [ones(1,size([sim_x sim_t1 sim_t2],2)) Inf];
        x0 = [.5*ones(1,size([sim_x sim_t1 sim_t2],2)) 1];
        [inp,~,~,~] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);
        EmulatorCovHypers = inp ;
    end
else
    if isequal(EmulatorMean,'Default')
        mean_sim = @(a,b,c) dual_calib_example_fn(...
            a,min_x,range_x,b,min_t1,range_t1,c,min_t2,range_t2,...
            mean_y,std_y);
    else
        mean_sim = EmulatorMean;
    end
    EmulatorCovHypers = [.5* ones(1,dim_x+dim_t1+dim_t2) Inf];
end

%% Set observation discrepancy mean
if isequal(obs_discrep_mean,'Default')
    mean_obs = @(x,t) zeros(size(x,1),1);
else
    mean_obs = obs_discrep_mean;
end

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
theta2_init = rand(dim_t2,1);
Sigma_theta1 = eye(size(theta1_init,1));
Sigma_theta2 = eye(size(theta2_init,1));


%% Set real and desired discrepancy rho and lambda prior distributions
log_obs_rho_prior  = @(r) sum(log( betapdf(r,2,0.4) ));
log_obs_lambda_prior = @(ld) log( gampdf(ld,150,4) );
log_des_rho_prior  = @(r) sum(log( betapdf(r,2,0.4) ));
log_des_lambda_prior = @(ld) log( gampdf(ld,150,4) );


%% Set real and desired discrepancy rho and lambda proposal distributions
% if isequal(obs_rho_lambda,'Default')
% We'll draw logit-transformed rho and lambda from normals centered at
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
if obs_discrep
    obs_rho_init = rand(dim_x+dim_t2,1);
    obs_lambda_init = gamrnd(1,1);
else
    obs_rho_init = .5*ones(dim_x+dim_t2,1);
    obs_lambda_init = Inf;
end
if des_discrep
    des_rho_init = rand(size(des_x,2),1);
    des_lambda_init = gamrnd(1,1);
else
    des_rho_init = .5*ones(size(des_x,2),1);
    des_lambda_init = Inf;
end
obs_Sigma_rho = eye(size(obs_rho_init,1));
obs_Sigma_lambda = 1;
des_Sigma_rho = eye(size(des_rho_init,1));
des_Sigma_lambda = 1;

%% Make observation and target error/tolerance covariance matrices
num_obs = size(obs_y(:),1) ; 
obs_cov_mat = obs_var * eye(num_obs) ;
num_tgt = size(des_y(:),1) ;
des_cov_mat = des_var * eye(num_tgt) ;

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
    'des_cov_mat',des_cov_mat,...
    'additional_discrep_cov',additional_discrep_cov,...
    'additional_discrep_mean',additional_discrep_mean,...
    'mean_sim',mean_sim,...
    'emulator_rho',EmulatorCovHypers(1:end-1),...
    'emulator_lambda',EmulatorCovHypers(end),...
    'mean_obs',mean_obs,...
    'rho_proposal',rho_proposal,...
    'lambda_proposal',lambda_proposal,...
    'rho_prop_log_mh_correction',rho_prop_log_mh_correction,...
    'lambda_prop_log_mh_correction',lambda_prop_log_mh_correction,...
    'des_rho_init',des_rho_init,...
    'des_lambda_init',des_lambda_init,...
    'des_rho_prop_cov',des_Sigma_rho,...
    'des_lambda_prop_cov',des_Sigma_lambda,...
    'obs_rho_init',obs_rho_init,...
    'obs_lambda_init',obs_lambda_init,...
    'obs_rho_prop_cov',obs_Sigma_rho,...
    'obs_lambda_prop_cov',obs_Sigma_lambda,...
    'log_obs_rho_prior',log_obs_rho_prior,...
    'log_obs_lambda_prior',log_obs_lambda_prior,...
    'log_des_rho_prior',log_des_rho_prior,...
    'log_des_lambda_prior',log_des_lambda_prior,...
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
    'obs_discrep',obs_discrep,...
    'des_discrep',des_discrep,...
    'modular',modular,...
    'doplot',doplot);

end