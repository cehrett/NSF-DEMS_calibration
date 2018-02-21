% Using samples on standardized scale:

post_mean = mean(samples(burn_in:end,:)) 

post_mean_input = [obs_x repmat(post_mean,size(obs_x,1),1) ];
em = emulator(tdat.input,tdat.output,post_mean_input,omega,rho,lambda,...
    0,0,true);
post_mean_sim_defl_std = em.mu(1:(length(em.mu)/3));
post_mean_sim_rot_std  = em.mu((length(em.mu)/3+1):(2*length(em.mu)/3));
post_mean_sim_cost_std = em.mu((2*length(em.mu)/3+1):end);
post_mean_sim_defl = post_mean_sim_defl_std * mean(tdat.output_sds(1,:))...
    + mean(tdat.output_means(1,:));
post_mean_sim_rot  = post_mean_sim_rot_std  * mean(tdat.output_sds(2,:))...
    + mean(tdat.output_means(2,:));
post_mean_sim_cost = post_mean_sim_cost_std * mean(tdat.output_sds(3,:))...
    + mean(tdat.output_means(3,:));
mean(post_mean_sim_defl) % range: .65-.83
mean(post_mean_sim_rot)  % range: 0.077-0.1
mean(post_mean_sim_cost) % range: 96-352

%% Using samples on original scale

vf_min = tdat.input_mins(3); 
vf_range = tdat.input_ranges(3);
thk_min = tdat.input_mins(4);
thk_range = tdat.input_ranges(4);

post_mean = [0.5879 16.7603];
post_mean = [0.5772 23.4230];
post_mean = [0.5577 22.6098];
post_mean = (post_mean - [vf_min thk_min] ) ./ [vf_range thk_range]

post_mean_input = [obs_x repmat(post_mean,size(obs_x,1),1) ];
em = emulator(tdat.input,tdat.output,post_mean_input,omega,rho,lambda,...
    0,0,true);
post_mean_sim_defl_std = em.mu(1:(length(em.mu)/3));
post_mean_sim_rot_std  = em.mu((length(em.mu)/3+1):(2*length(em.mu)/3));
post_mean_sim_cost_std = em.mu((2*length(em.mu)/3+1):end);
post_mean_sim_defl = post_mean_sim_defl_std * mean(tdat.output_sds(1,:))...
    + mean(tdat.output_means(1,:));
post_mean_sim_rot  = post_mean_sim_rot_std  * mean(tdat.output_sds(2,:))...
    + mean(tdat.output_means(2,:));
post_mean_sim_cost = post_mean_sim_cost_std * mean(tdat.output_sds(3,:))...
    + mean(tdat.output_means(3,:));
mean(post_mean_sim_defl) % range: .65-.83
mean(post_mean_sim_rot)  % range: 0.077-0.1
mean(post_mean_sim_cost) % range: 96-352

