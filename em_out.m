function output = em_out(samples,settings)
% This function takes as input the results of an MCMC routine and returns
% the output of the emulator at the posterior mean of the MCMC draws

% Unpack variables for convenience
burn_in = settings.burn_in;
obs_x = settings.obs_x;
tdat.input = settings.sim_xt;
tdat.output= settings.eta;
tdat.output_sds = settings.output_sds;
tdat.output_means = settings.output_means;
omega = settings.omega;
rho = settings.rho;
lambda = settings.lambda;
which_outputs = settings.which_outputs;

% Using samples on standardized scale

% How many outputs are we working with
num_out = sum(which_outputs);

post_mean = mean(samples(burn_in:end,:)) 

post_mean_input = [obs_x repmat(post_mean,size(obs_x,1),1) ];
em = emulator(tdat.input,tdat.output,post_mean_input,omega,rho,lambda,...
    0,0,true);

if num_out == 3 
    post_mean_sim_defl_std =em.mu(1:(length(em.mu)/3));
    post_mean_sim_rot_std  =em.mu((length(em.mu)/3+1):(2*length(em.mu)/3));
    post_mean_sim_cost_std =em.mu((2*length(em.mu)/3+1):end);
    post_mean_sim_defl = ...
        post_mean_sim_defl_std * mean(tdat.output_sds(1,:)) + ...
        mean(tdat.output_means(1,:));
    post_mean_sim_rot  = ...
        post_mean_sim_rot_std  * mean(tdat.output_sds(2,:)) + ...
        mean(tdat.output_means(2,:));
    post_mean_sim_cost = ...
        post_mean_sim_cost_std * mean(tdat.output_sds(3,:)) + ...
        mean(tdat.output_means(3,:));
    output = [ mean(post_mean_sim_defl)   ; % range: .65-.83
               mean(post_mean_sim_rot)    ; % range: 0.077-0.1
               mean(post_mean_sim_cost) ] ; % range: 96-352
else
    post_mean_sim_defl_std = em.mu(1:(length(em.mu)/2));
    post_mean_sim_cost_std = em.mu((length(em.mu)/2+1):end);
    
    post_mean_sim_defl = ...
        post_mean_sim_defl_std * mean(tdat.output_sds(1,:)) + ...
        mean(tdat.output_means(1,:));
    post_mean_sim_cost = ...
        post_mean_sim_cost_std * mean(tdat.output_sds(2,:))...
        + mean(tdat.output_means(2,:));
    output = [ mean(post_mean_sim_defl)   ; % range: .65-.83
               mean(post_mean_sim_cost) ] ; % range: 96-352
end

end