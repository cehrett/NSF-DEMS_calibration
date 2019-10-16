function settings = MCMC_settings (desired_obs,sim_x,sim_t,sim_y,varargin)
% settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...).
% 
% Parameter/Value pairs:
% 'M':
%   Number of samples to draw in the MCMC, includes burn-in. Default: 1e4.
% 'burn_in':
%   Proportion of samples to treat as burn-in. Default: 1/2.
% 'ObsVar':
%   'RefPrior' - (Default) Puts independent 1/sigma^2 prior on each obs var
%                          for each of the outputs.
%   'STOV'     - Uses set total observation variance. Default 10.
%   'Constant' - Uses constant observation variance. Default 0.05.
% 'ObsVarLvl':
%   Scalar value which gives the variance of the (constant) observation
%   error. Only used in conjunction with option ObsVar set to 'Constant'.
%   Default: 0.05.
% 'Cost_lambda':
%   false      - (Default) does not place informative prior on vf, thk.
%   true       - Places exp(-||(vf,thk)||^2) prior on vf, thk.
% 'which_outputs':
%   Row vector of 0/1 values indicating which of the outputs are to be
%   used. Default is to use all inputs.
% 'Rho_lam_optimum':
%   Row vector: [omega rho lambda]. Default is:
%               [   0.280981573480363   0.999189406633873...
%               0.600440750045477  0.719652153362981   0.102809702497319...
%               0.000837772517865 ]. If set to 0, then ML estimation via
%   gradient descent is used to find appropriate values.
% 'Discrepancy':
%   false      - (Default) No discrepancy function is used.
%   true       - Discrepancy function is used.
% 'DiscMargPrecProp':
%   Function of two inputs which serves as the proposal density for the
%   marginal precision on the discrepancy function. Default: identity
%   function on input 1. Also popular: @(x,s) exp(mvnrnd(log(x),s)).
% 'DiscMargPrecLogMHCorr':
%   Function of two inputs which serves as the log Metropolis-Hastings
%   correction factor for calculating the acceptance ratio in MCMC.
%   Default: constant zero function. If using the proposal function given
%   by @(x,s) exp(mvnrnd(log(x),s)), then use:
%   @(sig_s,sig)log(prod(sig_s))-log(prod(sig))
% 'LambdaDeltaInit':
%   Initial value for lambda_delta, the marginal precision for the
%   discrepancy function. Default: 1. Also popular: gamrnd(1,1).
% 'doplot':
%   true       - (Default) Update scatterplot every ten draws.
%   false      - No plots during MCMC.

%% Parse inputs
p = inputParser;
p.addRequired('desired_obs',@ismatrix);
p.addRequired('sim_x',@ismatrix);
p.addRequired('sim_t',@ismatrix);
p.addRequired('sim_y',@ismatrix);
p.addParameter('M',1e4,@isscalar);
p.addParameter('burn_in',1/2,@(x) x>0 && x<1);
p.addParameter('ObsVar','RefPrior',@isstr);
p.addParameter('Cost_lambda',0,@isscalar);
p.addParameter('which_outputs',ones(size(desired_obs)),@ismatrix);
p.addParameter('Rho_lam_optimum',...
        [0 999],... 
        @ismatrix);
p.addParameter('Discrepancy',false,@islogical);
p.addParameter('doplot',true,@islogical);
p.addParameter('DiscMargPrecProp',@(x,s)x,@(x)isa(x,'function_handle'));
p.addParameter(...
    'DiscMargPrecLogMHCorr',@(x,y)0,@(x)isa(x,'function_handle'));
p.addParameter('LambdaDeltaInit',1,@isscalar);
p.addParameter('ObsVarLvl',0.05,@isscalar);
p.parse(desired_obs,sim_x,sim_t,sim_y,varargin{:});

%% Collect inputs
M                    = p.Results.M;
burn_in              = floor(p.Results.burn_in*M);
ObsVar               = p.Results.ObsVar;
Cost_lambda          = p.Results.Cost_lambda;
which_outputs        = p.Results.which_outputs;
Rho_lam_optimum      = p.Results.Rho_lam_optimum;
Discrepancy          = p.Results.Discrepancy;
lambda_delta_init    = p.Results.LambdaDeltaInit;
lambda_prop_density  = p.Results.DiscMargPrecProp;
log_mh_correction_ld = p.Results.DiscMargPrecLogMHCorr;
ObsVarLvl            = p.Results.ObsVarLvl;
doplot               = p.Results.doplot;

%% Output reminder about Rho_lam_optimum
if isequal(Rho_lam_optimum,[0 999])
%     fprintf(['Using hyperparameter MLEs previously estimated for '...
%         'toy simulation problem.\n']);
%     Rho_lam_optimum = ...
%         [   0.280981573480363   0.999189406633873   0.600440750045477...
%         0.719652153362981   0.102809702497319   0.000837772517865 ] ; 
    fprintf(['Using hyperparameter MLEs previously estimated for '...
        'wind turbine application.\n']);
    Rho_lam_optimum  = ...
        [   0.935753521438069   0.650946653103927   0.673593619101900...
        0.479684392594821   0.967330479380613   0.015203646313917 ] ;
end

%% Infer useful values
num_cal = size(sim_t,2);
num_out = size(sim_y,2);
num_cntrl = size(sim_x,2) + num_out - 1 ; % The num_out - 1 is for dum vars
num_cntrl_wod = size(sim_x,2) ; % Without dummy vars


%% Proposal density
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x)./(1+exp(x));
prop_density = @(x,Sigma) logit_inv(mvnrnd(logit(x),Sigma)); % Normal
Sigma = [.5 0 ; 0 .5]; % Initial variance for prop_density
% Upper and lower bounds for theta (if applicable):
LB = min(sim_t) ; UB = max(sim_t) ; % Here it is assumed that the code was
                                    % run at extreme theta values
% Nugget size for computational stability of covariance matrices
nugsize = @(Covmat) 1e-4 ; % Why make it a function? So later it can be 
                           % made fancier if desired.
                           
%% MH correction for using log-normal proposal
log_mh_correction = @(theta_s,theta) log(prod(theta_s)*prod(1-theta_s))-...
    log(prod(theta)*prod(1-theta));
                           
%% Set prior for sigma^2_y: 1/sigma^2_y, and log version
switch ObsVar
    case 'RefPrior'
        % Joint multivariate draw version:
        sigma2=rand(1,num_out)*20;
        log_sigma2_prior = @(sigma2) -log(prod(sigma2));
        sigma2_prop_density = @(x,s) exp(mvnrnd(log(x),s));
        log_sig_mh_correction=@(sig_s,sig)log(prod(sig_s))-log(prod(sig));
        Sigma_sig = eye(num_out);
        init_sigma2_divs = 'null';
    case 'STOV'
        sigma2 = 10; % Sum total observation variance
        log_sigma2_prior = @(x) 0 ; % No prior on total variance
        init_sigma2_divs = [ 1/3 2/3 ] ; % Initial weights
        log_sig_mh_correction = @(x,s) 0 ; % No MH correction for sigma2
        sigma2_prop_density = @(x,s) x + rand(1,2) * s - s/2;
        Sigma_sig = .1; % Width of unif proposal window
        init_sigma2_divs = 0:1/num_out:1 ; % Initial proportion of total ov
        init_sigma2_divs = init_sigma2_divs(2:(end-1));
    case 'Constant'
        sigma2 = ObsVarLvl * ones(1,num_out); % Obs. var. of each output
        log_sigma2_prior = @(x) 0; % No prior on obs var.
        log_sig_mh_correction = @(x,s) 0; % No MH correction needed
        sigma2_prop_density = @(x,s) x;
        Sigma_sig = 0; % No var on proposal for sigma2
        init_sigma2_divs = 'null';
    otherwise
        error('InptErr','You specified an invalid value for ObsVar.');
end

%% Set prior for theta
Cost            = @(t,Cost_lambda) Cost_lambda * norm(t)^2;
theta_prior     = @(theta,Cost_lambda) exp(-Cost(theta,Cost_lambda));
log_theta_prior = @(theta,Cost_lambda) -Cost(theta,Cost_lambda);

%% Set initial discrepancy covariance parameters and priors
if Discrepancy 
    % Set initial values, proposal densities, prior densities, and any
    % requisite Metropolis Hastings corrections needed for omega_delta and
    % lambda_delta, which are the hyperparameters of the discrepancy
    % function.
    omega_delta_init       = betarnd(1,0.3,1,num_cntrl); % Initial value
    omega_prop_density     = @(x,Sigma) logit_inv(mvnrnd(logit(x),Sigma)); 
    Sigma_od               = .5*eye(num_cntrl);% Initial var for prop_dens
    log_mh_correction_od   = @(od_s,od) log(prod(od_s)*prod(1-od_s))-...
        log(prod(od)*prod(1-od));% log Metr. Hast. correction for proposal
    Sigma_ld               = 1 ; % Var for prop density on lambda_delta
    log_omega_delta_prior  = @(od) sum(log( betapdf(od,1,0.6) ));
    log_lambda_delta_prior = @(ld) log( gampdf(ld,5,5) );
else 
    % If we're not using a discrepancy function, then set all these to
    % accordingly.
    omega_delta_init       = .5 * ones(1,num_cntrl) ;
    lambda_delta_init      = Inf;
    omega_prop_density     = @(x,Sigma)x;
    Sigma_od               = 'null';
    log_mh_correction_od   = 'null';
    lambda_prop_density    = @(x)x;
    Sigma_ld               = 'null';
    log_mh_correction_ld   = 'null';
    log_omega_delta_prior  = @(od)0;
    log_lambda_delta_prior = @(ld)0;
end

%% Package proposal density
% Here, collect the proposal densities and hyperparamters for theta, for
% omega_delta, and for lambda_delta, along with any requisite
% Metropolis-Hastings corrections.
proposal.density              = prop_density; 
proposal.Sigma                = Sigma;
proposal.sigma2_prop_density  = sigma2_prop_density;
proposal.Sigma_sig            = Sigma_sig;
proposal.omega_prop_density   = omega_prop_density;
proposal.lambda_prop_density  = lambda_prop_density;
proposal.Sigma_od             = Sigma_od;
proposal.log_mh_correction_od = log_mh_correction_od;
proposal.Sigma_ld             = Sigma_ld;
proposal.log_mh_correction_ld = log_mh_correction_ld;


%% Load data and get initial theta value
raw_dat = [sim_x sim_t sim_y];
indx = 1:num_out; % This will help tell which columns of raw_dat we need.
for ii = 1:length(which_outputs) % This loop will set indx appropriately.
    if which_outputs(ii) indx = [indx num_cal+num_cntrl_wod+ii] ; end
end
raw_dat = raw_dat(:,indx);

% Get initial theta val
init_theta = rand(1,num_cal) 

%% Rescale inputs, standardize outputs
tdat = Tdat(raw_dat,num_out); % rescaling inputs, standardizing outputs. 
                              % Tdat also adds dummy control variables
                              % with num_out levels.
fprintf('done.\n\n')
sim_xt = tdat.input; % Normalized model inputs
eta = tdat.output;   % Standardized model outputs
% Input mins, input ranges, output means, and output sd's are stored so
% that the rescaling and standardization can be reversed later.
sim_cntrl_input_mins = tdat.input_mins(1:num_cntrl) ;
sim_cntrl_input_ranges = tdat.input_ranges(1:num_cntrl) ;
sim_calib_input_mins = tdat.input_mins(num_cntrl+1:end) ;
sim_calib_input_ranges = tdat.input_ranges(num_cntrl+1:end) ;
sim_output_means = mean(tdat.output_means') ;
% Note that in the above line, the mean is taken across temperatures. This
% is because the standardization was carried out by individual temperature
% setting, whereas since our desired data will be constant across
% temperature, we want a single value for all temperature settings. This
% fudge should be harmless, in that the values do not differ much across
% temperature. Sd's will be treated similarly in the following line.
sim_output_sds = mean(tdat.output_sds') ;
% Assuming uniform prior for theta: find LB,UB and rescale them
LB = (LB - sim_calib_input_mins)./sim_calib_input_ranges ;
UB = (UB - sim_calib_input_mins)./sim_calib_input_ranges ;

%% MCMC settings
% Covariance parameter settings, found by optimization routine:
if Rho_lam_optimum == 0
    fprintf('No rho,lambda values specified; commencing ML estimation.\n');
    Rho_lam_optimum = opt_rho_lambda(sim_xt,eta,num_cal,...
        rand(1,num_cntrl),rand(1,num_cal),gamrnd(5,5));
end
omega  = Rho_lam_optimum(1:num_cntrl);
rho    = Rho_lam_optimum(num_cntrl+1:num_cntrl+num_cal);
lambda = Rho_lam_optimum(num_cntrl+num_cal+1);
% Need different omega in case not all three outputs used
if sum(which_outputs) ~= 3
    error('To use < 3 outputs, sort out manually what to do with omega');
end

%% Prepare field observations, control settings, and variance for MCMC
% First, we must standardize them to be on the same scale as the simulation
% observations.
y = (desired_obs - sim_output_means) ./ sim_output_sds ;
% Now, pair these standardized observations with appropriate control
% settings. Notice that we are here assuming only one field observation,
% which is constant across control settings (other than dummy variable).
obs_x = unique(tdat.input(:,1:num_cntrl),'rows','stable');
y = repelem(y,size(obs_x,1)/length(y))' ;


%% Set uniform prior for theta
out_of_range = @(theta) theta < LB | theta > UB ; 

settings = struct(...
    'M',M,...
    'burn_in',burn_in,...
    'sim_xt',sim_xt,...
    'num_cntrl',num_cntrl,...
    'num_cntrl_wod',num_cntrl_wod,...
    'num_cal',num_cal,...
    'num_out',num_out,...
    'eta',eta,...
    'obs_x',obs_x,...
    'y',y,...
    'sigma2',sigma2,...
    'init_sigma2_divs',init_sigma2_divs,...
    'log_sigma2_prior',log_sigma2_prior,...
    'log_omega_delta_prior',log_omega_delta_prior,...
    'log_lambda_delta_prior',log_lambda_delta_prior,...
    'omega_delta_init',omega_delta_init,...
    'lambda_delta_init',lambda_delta_init,...
    'out_of_range',out_of_range,...
    'init_theta',init_theta,...
    'omega',omega,...
    'rho',rho,...
    'lambda',lambda,...
    'proposal',proposal,...
    'nugsize',nugsize,...
    'log_sig_mh_correction',log_sig_mh_correction,...
    'log_mh_correction',log_mh_correction,...
    'input_cntrl_mins',min(sim_x),...
    'input_calib_mins',min(sim_t),...
    'input_cntrl_ranges',range(sim_x),...
    'input_calib_ranges',range(sim_t),...
    'output_sds',tdat.output_sds,...
    'output_means',tdat.output_means,...
    'log_theta_prior',log_theta_prior,...
    'Cost_lambda',Cost_lambda,...
    'which_outputs',which_outputs,...
    'desired_obs',desired_obs,...
    'doplot',true);

end