function emout = em_out_many(samples,settings,n)
% This function takes as input the results of an MCMC routine and returns
% the output of the emulator at the locations drawn in the MCMC.
% This function is a modified version of em_out. em_out is designed to find
% the output of the posterior mean from the MCMC. By constrast, em_out_many
% is designed to give the model output at many or all draws from the MCMC.
% The parameter n gives the desired number of samples at which to find the
% output. n instances will be taken from samples. To use all (post-burn in)
% samples, set n=0.

% Unpack variables for convenience
burn_in = settings.burn_in;
tdat.input = settings.sim_xt;
tdat.output= settings.eta;
tdat.output_sds = settings.output_sds;
tdat.output_means = settings.output_means;
omega = settings.omega;
rho = settings.rho;
lambda = settings.lambda;
which_outputs = settings.which_outputs;

% Following line reduces the set of samples down to a size-n set of
% post-burn-in samples, if n!=0; otherwise all post-burn_in samples used
samples = samples(burn_in:end,:);
if n ~= 0
    samples = datasample(samples,n,1,'Replace',false);
end

% Get control settings
dum_vars = unique(settings.obs_x(:,1));
obs_x = repmat(dum_vars,size(samples,1),1);
obs_x = [ obs_x 0.5 * ones(size(obs_x,1),1) ];

% How many outputs are we working with
num_out = sum(which_outputs);

% Get prediction locations
pred_locs = [obs_x repelem(samples,length(dum_vars),1) ];
em = emulator(tdat.input,tdat.output,pred_locs,omega,rho,lambda,...
    0,0,true);

% Get outputs and related stats. Note that the code here assumes that we
% are getting all three outputs of deflection, rotation, cost.
outputs_defl_mean_std = em.mu(obs_x(:,1)==0);
outputs_rotn_mean_std = em.mu(obs_x(:,1)==0.5);
outputs_cost_mean_std = em.mu(obs_x(:,1)==1);
outputs_defl_sd_std = sqrt(diag(em.cov_gp(obs_x(:,1)==0,obs_x(:,1)==0)));
outputs_rotn_sd_std = sqrt(diag(em.cov_gp(obs_x(:,1)==.5,obs_x(:,1)==.5)));
outputs_cost_sd_std = sqrt(diag(em.cov_gp(obs_x(:,1)==1,obs_x(:,1)==1)));

% Transform outputs back to original scale
outputs_defl_mean =outputs_defl_mean_std * mean(tdat.output_sds(1,:)) + ...
    mean(tdat.output_means(1,:));
outputs_rotn_mean =outputs_rotn_mean_std * mean(tdat.output_sds(2,:)) + ...
    mean(tdat.output_means(2,:));
outputs_cost_mean =outputs_cost_mean_std * mean(tdat.output_sds(3,:)) + ...
    mean(tdat.output_means(3,:));
outputs_defl_sd = outputs_defl_sd_std * mean(tdat.output_sds(1,:)) + ...
    mean(tdat.output_means(1,:));
outputs_rotn_sd = outputs_rotn_sd_std * mean(tdat.output_sds(2,:)) + ...
    mean(tdat.output_means(2,:));
outputs_cost_sd = outputs_cost_sd_std * mean(tdat.output_sds(3,:)) + ...
    mean(tdat.output_means(3,:));
output_means = [outputs_defl_mean outputs_rotn_mean outputs_cost_mean];
output_sds   = [outputs_defl_sd outputs_rotn_sd outputs_cost_sd];
emout = struct('output_means',output_means,'output_sds',output_sds);
% range: .65-.83
% range: 0.077-0.1
% range: 96-352

end