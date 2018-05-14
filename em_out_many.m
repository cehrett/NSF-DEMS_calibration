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
num_cntrl = settings.num_cntrl;

% Following line reduces the set of samples down to a size-n set of
% post-burn-in samples, if n!=0; otherwise all post-burn_in samples used
samples = samples(burn_in:end,:);
if n ~= 0
    samples = datasample(samples,n,1,'Replace',false);
end

% Get control settings
num_cntrl_wodv = num_cntrl - ( num_out - 1 ) ; % # cntrl vars w/o dum vars
dum_vars = unique(settings.obs_x(:,1:end-num_cntrl_wodv),'stable');
obs_x = repmat(dum_vars,size(samples,1),1);
% Here is where the observation points are put together. The final column
% of obs_x corresponds to the non-dummy variable control variables. In the
% current version of the code, they are simply set to be equal to the
% midpoint of the range of each such variable. If these control variables
% are fairly important, then this will be a bad idea, since it will ignore
% the variability in this/these control variable/s. 
obs_x = [ obs_x 0.5 * ones(size(obs_x,1),num_ctrl_wodv) ];

% How many outputs are we working with
num_out = sum(which_outputs);

% Get prediction locations
pred_locs = [obs_x repelem(samples,length(dum_vars),1) ];
em = emulator(tdat.input,tdat.output,pred_locs,omega,rho,lambda,...
    0,0,true);




% Get outputs and transform outputs back to original scale
% initialize variables
outputs_defl_mean = [];
outputs_rotn_mean = [];
outputs_cost_mean = [];
outputs_defl_sd   = [];
outputs_rotn_sd   = [];
outputs_cost_sd   = [];
ndx = 1; % Just to tell us which output we're working with
if which_outputs(1)==1
    ndcs = obs_x(:,ndx) == 1; % Currently relevant indices
    outputs_defl_mean_std = em.mu(ndcs);
    outputs_defl_sd_std=sqrt(diag(em.cov_gp(ndcs,ndcs)));
    outputs_defl_mean =outputs_defl_mean_std * ...
        mean(tdat.output_sds(ndx,:)) + ...
        mean(tdat.output_means(ndx,:));
    outputs_defl_sd = outputs_defl_sd_std * mean(tdat.output_sds(ndx,:));
    ndx = ndx + 1;
end
if which_outputs(2)==1
    ndcs = obs_x(:,ndx) == 1; % Currently relevant indices
    outputs_rotn_mean_std = em.mu(ndcs);
    outputs_rotn_sd_std=sqrt(diag(em.cov_gp(ndcs,ndcs)));
    outputs_rotn_mean =outputs_rotn_mean_std * ...
        mean(tdat.output_sds(ndx,:)) + ...
        mean(tdat.output_means(ndx,:));
    outputs_rotn_sd = outputs_rotn_sd_std * mean(tdat.output_sds(ndx,:));
    ndx = ndx + 1;
end
if which_outputs(3)==1
    ndcs = all(obs_x(:,1:ndx-1)==zeros(1,ndx-1),2); % Currently rel indices
    outputs_cost_mean_std = em.mu(ndcs);
    outputs_cost_sd_std = sqrt(diag(em.cov_gp(ndcs,ndcs)));
    outputs_cost_mean =outputs_cost_mean_std * ...
        mean(tdat.output_sds(ndx,:)) + ...
        mean(tdat.output_means(ndx,:));
    outputs_cost_sd = outputs_cost_sd_std * mean(tdat.output_sds(ndx,:)); 
end
 
 

output_means = [outputs_defl_mean outputs_rotn_mean outputs_cost_mean];
output_sds   = [outputs_defl_sd outputs_rotn_sd outputs_cost_sd];

% Pack up and leave
emout = struct('output_means',output_means,'output_sds',output_sds);


end