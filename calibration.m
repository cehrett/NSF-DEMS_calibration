% Master file

%% Add paths
addpath('.\NSF DEMS\Phase 1\');
addpath('.\NSF DEMS\Phase 1\stored_data');

%% USER SET VALUES
% Set desired data:
desired_obs = [.65, .077, 96] ; % tip deflection, rotation, and cost
desired_sds = desired_obs / 4 ; 
% MCMC settings
M = 1e3 ; % Total number of draws (including burn-in)
burn_in = ceil(M/5) ; % burn-in
% Proposal density
prop_density = @(x,Sigma) (mvnrnd(x,Sigma)); % Normal
Sigma = [.05 0 ; 0 .05]; % Initial variance for prop_density
% Upper and lower bounds for theta (if applicable)
LB = [ .2 10 ] ; UB = [.6 25] ;
% Nugget size for computational stability of covariance matrices
nugsize = @(Covmat) 1e-4 ; % Why make it a function? So later it can be 
                           % made fancier.
% Covariance parameter settings, found by optimization routine:
Rho_lam_optimum  = [0.655344235568109   0.931941001705886 ...
    0.960653924901867   0.991953924787049   0.017385994893994];
omega  = Rho_lam_optimum(1:2);
rho    = Rho_lam_optimum(3:4);
lambda = Rho_lam_optimum(5);
%% END USER SET VALUES

%% Package proposal density
proposal.density = prop_density; proposal.Sigma = Sigma;

%% Load data and get initial theta value
fprintf('Reading data from .xlsx...\n')
raw_dat = xlsread("fe_results.xlsx");
% raw_dat is assumed to include one observation per row, with input columns
% preceding output columns (with no headers).
num_out = 3; % This is the number of outputs in the data, which are assumed
             % to be the final columns of raw_dat.
num_calib = 2 ; % This is the number of calibration parameters in the data,
                % which are assumed to be the columns immediately preceding
                % the output columns in raw_dat.
num_cntrl = size(raw_dat,2) - num_out - num_calib + (num_out>1) ; 
% num_cntrl is the number of control inputs in raw_dat, plus one if the
% output is multivariate; the extra one is for a dummy input.
init_theta = rand(1,2) 

%% Rescale inputs, standardize outputs
tdat = Tdat(raw_dat,num_out); % rescaling inputs, standardizing outputs. 
                              % Tdat also adds a dummy control variable
                              % with num_out levels.
fprintf('done.\n\n')
sim_xt = tdat.input ;
eta = tdat.output;
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

%% Prepare field observations, control settings, and variance for MCMC
% First, we must standardize them to be on the same scale as the simulation
% observations.
y = (desired_obs - sim_output_means) ./ sim_output_sds ;
% Now, pair these standardized observations with appropriate control
% settings. Notice that we are here assuming only one field observation,
% which is constant across control settings (other than dummy variable).
obs_x = unique(tdat.input(:,1:num_cntrl),'rows');
y = repelem(y,size(obs_x,1)/length(y))' ;
% Now set the observation variance matrix. First, standardize the sd's
% which are specified by the user.
y_sds = desired_sds ./ sim_output_sds ;
% Now make it into a covariance matrix:
y_sds = repelem(y_sds,size(obs_x,1)/length(y_sds))';
Sigma_y = diag(y_sds.^2);

%% Set uniform prior for theta
out_of_range = @(theta) theta < LB | theta > UB ; 

%% Run MCMC routine
[samples,Sigma] = MCMC_joint_proposal(M,burn_in,sim_xt,eta,obs_x,y,...
    Sigma_y,out_of_range,init_theta,omega,rho,lambda,proposal,nugsize)


