function output = em_out(samples,burn_in,obs_x,tdat_input,tdat_output,...
    tdat_output_sds,tdat_output_means,omega,rho,lambda,which_outputs)
% This function takes as input the results of an MCMC routine and returns
% the output of the emulator at the posterior mean of the MCMC draws

% tdat.input = tdat_input;
% tdat.output= tdat_output;
% tdat.output_sds = tdat_output_sds;
% tdat.output_means = tdat_output_means;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_out=3;
raw_dat = xlsread('fe_results.xlsx');
% raw_dat is assumed to include one observation per row, with input columns
% preceding output columns (with no headers).
num_calib = 2 ; % This is the number of calibration parameters in the data,
                % which are assumed to be the columns immediately preceding
                % the output columns in raw_dat.
num_cntrl = size(raw_dat,2) - num_out - num_calib + (num_out>1) ; 
% num_cntrl is the number of control inputs in raw_dat, plus one if the
% output is multivariate; the extra one is for a dummy input.

%% Rescale inputs, standardize outputs
tdat = Tdat(raw_dat,num_out); % rescaling inputs, standardizing outputs. 
                              % Tdat also adds a dummy control variable
                              % with num_out levels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using samples on standardized scale

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