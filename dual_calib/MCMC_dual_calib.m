function results = MCMC_dual_calib(settings)
% results = MCMC_dual_calib(settings)
%
% settings is a structure generated by the file MCMC_dual_calib_settings.m.


%% Unpack settings
M                             = settings.M;
burn_in                       = settings.burn_in;
sim_x                         = settings.sim_x_01;
sim_t1                        = settings.sim_t1_01;
sim_t2                        = settings.sim_t2_01;
sim_y                         = settings.sim_y_std;
obs_x                         = settings.obs_x_01;
obs_t2                        = settings.obs_t2_01;
obs_y                         = settings.obs_y_std;
des_x                         = settings.des_obs_x_01;
des_y                         = settings.des_obs_y_std;
min_x                         = settings.min_x;
range_x                       = settings.range_x;
min_t1                        = settings.min_t1;
range_t1                      = settings.range_t1;
min_t2                        = settings.min_t2;
range_t2                      = settings.range_t2;
mean_y                        = settings.mean_y;
std_y                         = settings.std_y;
obs_cov_mat                   = settings.obs_cov_mat;
emulator_rho                  = settings.emulator_rho;
emulator_lambda               = settings.emulator_lambda;
obs_discrep_rho               = settings.obs_discrep_rho;
obs_discrep_lambda            = settings.obs_discrep_lambda;
rho                           = settings.des_discrep_rho_init;
lambda                        = settings.des_discrep_lambda_init;
rho_proposal                  = settings.rho_proposal;
lambda_proposal               = settings.lambda_proposal;
rho_prop_log_mh_correction    = settings.rho_prop_log_mh_correction;
lambda_prop_log_mh_correction = settings.lambda_prop_log_mh_correction;
Sigma_rho                     = settings.rho_prop_cov;
Sigma_lambda                  = settings.lambda_prop_cov;
log_rho_prior_fn              = settings.des_discrep_log_rho_prior;
log_lambda_prior_fn           = settings.des_discrep_log_lambda_prior;
theta1                        = settings.theta1_init;
theta2                        = settings.theta2_init;
theta1_proposal               = settings.theta1_proposal;
theta2_proposal               = settings.theta2_proposal;
theta1_prop_log_mh_correction = settings.theta1_prop_log_mh_correction;
theta2_prop_log_mh_correction = settings.theta2_prop_log_mh_correction;
Sigma_theta1                  = settings.theta1_prop_cov;
Sigma_theta2                  = settings.theta2_prop_cov;
log_theta1_prior_fn           = settings.log_theta1_prior;
log_theta2_prior_fn           = settings.log_theta2_prior;
doplot                        = settings.doplot;


%% Set some useful variables
% nugsize tells us what size nugget to add to a matrix for computational
% stability. Here, it's just a constant, but is introduced as a function so
% that it can easily be upgraded to something fancier if desired.
nugsize = @(X) 1e-4;
% D is the vector of all outputs from all sources.
D = [ sim_y ; obs_y ; des_y ] ;
opt_acc_rate = 0.234 ; % Acceptance rate treated as optimal in MCMC
upd = 100 ; % Tune adaptive covariances (and update plots) every upd loops
% The next four variables will be used when the values being sampled are
% not univariate, so that the plots cycle through which element of each
% vector is displayed.
col_theta1 = 0;
col_theta2 = 0;
col_rho    = 0;
col_lambda = 0;


%% Set mean functions for GPs
% Currently, constant mean zero is used for all GPs. However, that is
% implemented here in function form so that it can be changed easily in the
% future if desired.
mean_sim = @(a,b,c) zeros(size(a,1)); % Emulator mean
mean_obs = @(a,b) zeros(size(a,1)); % Observation discrepancy mean
mean_des = @(a) zeros(size(a,1)); 


%% Prepare that part of covariance matrix which will not change in MCMC
% The cov matrix we need, Sigma_z, can be decomposed so that this big part
% of it remains unchanged, so that we need calculate that only this once.
% Massive computation savings over getting Sigma_z from scratch each time:
Sigma_emulator_simsim = ...
    gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],[sim_x sim_t1 sim_t2],...
    emulator_lambda,false);
% Get computationally nice versions of it and its inverse for posterior
% predictions:
Sess = Sigma_emulator_simsim + ...
    eye(size(Sigma_emulator_simsim)) * nugsize(Sigma_emulator_simsim);
iSess = inv(Sess);


%% Initialize some variables for later use
startplot          = 10                           ;
accepted_theta1    = 0                            ; % accepted theta1
accepted_theta2    = 0                            ; % accepted theta2
accepted_rho       = 0                            ; % accepted rho
accepted_lambda    = 0                            ; % accepted lambda
msg                = 0                            ; % For console output
theta1_rec         = zeros(M,numel(theta1))       ;
theta2_rec         = zeros(M,numel(theta2))       ;
theta1_rec(1,:)    = theta1                       ;
theta2_rec(1,:)    = theta2                       ;
rho_rec            = zeros(M,numel(rho))          ;
rho_rec(1,:)       = rho                          ;
lambda_rec         = zeros(M,numel(lambda))       ;
lambda_rec(1,:)    = lambda                       ;
mult_theta1        = 10                           ; % mult. for proposal
mult_theta2        = 10                           ; 
mult_rho           = 10                           ; 
mult_lambda        = 10                           ; 


%% Get initial log likelihoods
% First, get the rest of the covariance matrix. We need to get the
% covariance matrices for the emulator of each combination of simulator
% runs, real observations, and desired observations; then we need to get
% the covariance matrices for the true discrepancy for each combination of
% true observations and desired observations; then we need to get the
% covariance matrix for the target discrepancy on the desired observations.
% Then we need to combine all these into one big covariance matrix.
Sigma_emulator_simobs = ...
    gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        emulator_lambda,false);
Sigma_emulator_simdes = ...
    gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
Sigma_emulator_obsobs = ...
    gp_cov(emulator_rho,...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        emulator_lambda,false);
Sigma_emulator_obsdes = ...
    gp_cov(emulator_rho,...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
Sigma_emulator_desdes = ...
    gp_cov(emulator_rho,...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
Sigma_obs_dscr_obsobs = ...
    gp_cov(obs_discrep_rho,[obs_x obs_t2],[obs_x obs_t2],...
        obs_discrep_lambda,false);
Sigma_obs_dscr_obsdes = ...
    gp_cov(obs_discrep_rho,[obs_x obs_t2],...
        [des_x repmat(theta2,size(des_x,1))],...
        obs_discrep_lambda,false);
Sigma_obs_dscr_desdes = ...
    gp_cov(obs_discrep_rho,[des_x repmat(theta2,size(des_x,1))],...
        [des_x repmat(theta2,size(des_x,1))],...
        obs_discrep_lambda,false);
Sigma_des_dscr_desdes = ...
    gp_cov(rho,des_x,des_x,lambda,false);
% Now to combine everything into one big covariance matrix:
Sigma = ...
    [Sigma_emulator_simsim Sigma_emulator_simobs' Sigma_emulator_simdes' ;
    Sigma_emulator_simobs ...
        Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + obs_cov_mat ...
        Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes';
    Sigma_emulator_simdes ...
        Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes ...
        Sigma_emulator_desdes+Sigma_obs_dscr_desdes+Sigma_des_dscr_desdes];
% Add nugget for computational tractability:
Sigma = Sigma + eye(size(Sigma)) * nugsize(Sigma) ;
% Now to get the log likelihoods
mu = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
    mean_sim(obs_x,repmat(theta1,size(obs_x,1)),obs_t2) + ...
        mean_obs(obs_x,obs_t2);
    mean_sim(des_x,repmat(theta1,size(des_x,1)),...
        repmat(theta2,size(des_x,1))) + ...
        mean_obs(des_x,repmat(theta2,size(des_x,1))) + ...
        mean_des(des_x)]; % Mean of D
log_cond_dens_D = logmvnpdf(D',mu,Sigma);
log_theta1_prior = log_theta1_prior_fn(theta1);
log_theta2_prior = log_theta2_prior_fn(theta2);
log_rho_prior = log_rho_prior_fn(rho);
log_lambda_prior = log_lambda_prior_fn(lambda);

%%%%%%%%%%%%%%%
%% MCMC loop %%
%%%%%%%%%%%%%%%

for ii = 2:M
    
    %% Draw new theta1
    theta1_s = theta1_proposal(theta1,Sigma_theta1);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw theta1_s, using the updated mean and
    % covariance of the GPs. Only certain parts of the big covariance
    % matrix Sigma need to be updated -- those that depend on theta1.
    mu_s = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
        mean_sim(obs_x,repmat(theta1_s,size(obs_x,1)),obs_t2) + ...
            mean_obs(obs_x,obs_t2);
        mean_sim(des_x,repmat(theta1_s,size(des_x,1)),...
            repmat(theta2,size(des_x,1))) + ...
            mean_obs(des_x,repmat(theta2,size(des_x,1))) + ...
            mean_des(des_x)]; % Mean of D
    Sigma_emulator_simobs_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],emulator_lambda,...
        false);
    Sigma_emulator_simdes_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [des_x repmat(theta1_s,size(obs_x,1),1) ...
            repmat(theta2,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_obsobs_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        emulator_lambda,false);
    Sigma_emulator_obsdes_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        [des_x repmat(theta1_s,size(obs_x,1),1) ...
            repmat(theta2,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_desdes_s = gp_cov(emulator_rho,...
        [des_x repmat(theta1_s,size(obs_x,1),1) ...
            repmat(theta2,size(obs_x,1),1)],...
        [des_x repmat(theta1_s,size(obs_x,1),1) ...
            repmat(theta2,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs_s' ...
            Sigma_emulator_simdes_s' ;
        Sigma_emulator_simobs_s ...
            Sigma_emulator_obsobs_s + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes_s' + Sigma_obs_dscr_obsdes';
        Sigma_emulator_simdes_s ...
            Sigma_emulator_obsdes_s + Sigma_obs_dscr_obsdes ...
            Sigma_emulator_desdes_s + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes...
    ];
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu_s,Sigma_s);
    log_theta1_prior_s = log_theta1_prior_fn(theta1_s);
    log_lik_theta1_s = log_cond_dens_D_s + log_theta1_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_theta1 = log_cond_dens_D + log_theta1_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_theta1_s - log_lik_theta1 + ...
        theta1_prop_log_mh_correction(theta1_s,theta1);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        theta1 = theta1_s;
        mu = mu_s;
        Sigma_emulator_simobs = Sigma_emulator_simobs_s;
        Sigma_emulator_simdes = Sigma_emulator_simdes_s;
        Sigma_emulator_obsobs = Sigma_emulator_obsobs_s;
        Sigma_emulator_obsdes = Sigma_emulator_obsdes_s;
        Sigma_emulator_desdes = Sigma_emulator_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_theta1 = accepted_theta1 + 1;
    end
        
    
    %% Draw new theta2
    theta2_s = theta2_proposal(theta2,Sigma_theta2);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw theta2_s, using the updated mean and
    % covariance of the GPs. Only certain parts of the big covariance
    % matrix Sigma need to be updated -- those that depend on theta2.
    mu_s = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
        mean_sim(obs_x,repmat(theta1,size(obs_x,1)),obs_t2) + ...
            mean_obs(obs_x,obs_t2);
        mean_sim(des_x,repmat(theta1,size(des_x,1)),...
            repmat(theta2_s,size(des_x,1))) + ...
            mean_obs(des_x,repmat(theta2_s,size(des_x,1))) + ...
            mean_des(des_x)]; % Mean of D
    Sigma_emulator_simdes_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [des_x repmat(theta1,size(obs_x,1),1) ...
            repmat(theta2_s,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_obsdes_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        [des_x repmat(theta1,size(obs_x,1),1) ...
            repmat(theta2_s,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_desdes_s = gp_cov(emulator_rho,...
        [des_x repmat(theta1,size(obs_x,1),1) ...
            repmat(theta2_s,size(obs_x,1),1)],...
        [des_x repmat(theta1,size(obs_x,1),1) ...
            repmat(theta2_s,size(obs_x,1),1)],...
        emulator_lambda,false);
    Sigma_obs_dscr_obsdes_s = gp_cov(obs_discrep_rho,...
        [obs_x obs_t2],[des_x repmat(theta2_s,size(obs_x,1),1)],...
        obs_discrep_lambda,false);
    Sigma_obs_dscr_desdes_s = gp_cov(obs_discrep_rho,...
        [des_x repmat(theta2_s,size(obs_x,1),1)],...
        [des_x repmat(theta2_s,size(obs_x,1),1)],...
        obs_discrep_lambda,false);
    Sigma_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs' ...
            Sigma_emulator_simdes_s' ;
        Sigma_emulator_simobs ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes_s' + Sigma_obs_dscr_obsdes_s';
        Sigma_emulator_simdes_s ...
            Sigma_emulator_obsdes_s + Sigma_obs_dscr_obsdes_s ...
            Sigma_emulator_desdes_s + Sigma_obs_dscr_desdes_s + ...
                Sigma_des_dscr_desdes...
    ];
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu_s,Sigma_s);
    log_theta2_prior_s = log_theta2_prior_fn(theta2_s);
    log_lik_theta2_s = log_cond_dens_D_s + log_theta2_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_theta2 = log_cond_dens_D + log_theta2_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_theta2_s - log_lik_theta2 + ...
        theta2_prop_log_mh_correction(theta2_s,theta2);
    
    % Now accept theta2_s with probability min(alpha,1)
    if log(rand) < log_alpha
        theta2 = theta2_s;
        mu = mu_s;
        Sigma_emulator_simdes = Sigma_emulator_simdes_s;
        Sigma_emulator_obsdes = Sigma_emulator_obsdes_s;
        Sigma_emulator_desdes = Sigma_emulator_desdes_s;
        Sigma_obs_dscr_obsdes = Sigma_obs_dscr_obsdes_s;
        Sigma_obs_dscr_desdes = Sigma_obs_dscr_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_theta2 = accepted_theta2 + 1;
    end
    
    %% Draw new rho
    rho_s = rho_proposal(rho,Sigma_rho);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw rho_s, using the updated covariance of
    % the GPs. Only part of the big covariance
    % matrix Sigma need to be updated -- the part that depends on rho.
    Sigma_des_dscr_desdes_s = gp_cov(rho_s,des_x,des_x,lambda,false);
    Sigma_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs' ...
            Sigma_emulator_simdes' ;
        Sigma_emulator_simobs ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes';
        Sigma_emulator_simdes ...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes_s...
    ];
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu,Sigma_s);
    log_rho_prior_s = log_rho_prior_fn(rho_s);
    log_lik_rho_s = log_cond_dens_D_s + log_rho_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_rho = log_cond_dens_D + log_rho_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_rho_s - log_lik_rho + ...
        rho_prop_log_mh_correction(rho_s,rho);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        rho = rho_s;
        Sigma_des_dscr_desdes = Sigma_des_dscr_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_rho = accepted_rho + 1;
    end
    
    %% Draw new lambda
    lambda_s = lambda_proposal(lambda,Sigma_lambda);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw lambda_s, using the updated covariance 
    % of the GPs. Only part of the big covariance
    % matrix Sigma need to be updated -- the part that depends on lambda.
    Sigma_des_dscr_desdes_s = gp_cov(rho,des_x,des_x,lambda_s,false);
    Sigma_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs' ...
            Sigma_emulator_simdes' ;
        Sigma_emulator_simobs ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes';
        Sigma_emulator_simdes ...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes_s...
    ];
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu,Sigma_s);
    log_lambda_prior_s = log_lambda_prior_fn(lambda_s);
    log_lik_lambda_s = log_cond_dens_D_s + log_lambda_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_lambda = log_cond_dens_D + log_lambda_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_lambda_s - log_lik_lambda + ...
        lambda_prop_log_mh_correction(lambda_s,lambda);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        lambda = lambda_s;
        Sigma_des_dscr_desdes = Sigma_des_dscr_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_lambda = accepted_lambda + 1;
    end
    
    %% Record draws
    theta1_rec(ii,:) = theta1;
    theta2_rec(ii,:) = theta2;
    rho_rec(ii,:)    = rho;
    lambda_rec(ii,:) = lambda;
    
    %% Adjust adaptive proposal covariances
    % We adjust the proposal covariances every upd draws.
    if mod(ii,upd) == 0 && ii <= burn_in
        
        fprintf(repmat('\b',1,msg)); % Prepare to output to console
        
        % Adjust proposal covariance for theta1
        mult_mult = max(.5,min(2,accepted_theta1/100/opt_acc_rate));
        mult_theta1 = mult_mult * mult_theta1;
        fprintf('theta1 proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_theta1 = 0 ;
        Sigma_theta1 = mult_theta1 * cov(logit(theta1_rec(1:ii,:)));
        
        % Adjust proposal covariance for theta2
        mult_mult = max(.5,min(2,accepted_theta2/100/opt_acc_rate));
        mult_theta2 = mult_mult * mult_theta2;
        fprintf('theta2 proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_theta2 = 0 ;
        Sigma_theta2 = mult_theta2 * cov(logit(theta2_rec(1:ii,:)));
        
        % Adjust proposal covariance for rho
        mult_mult = max(.5,min(2,accepted_rho/100/opt_acc_rate));
        mult_rho = mult_mult * mult_rho;
        fprintf('rho proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_rho = 0 ;
        Sigma_rho = mult_rho * cov(log(rho_rec(1:ii,:)));

        % Adjust proposal covariance for lambda
        mult_mult = max(.5,min(2,accepted_lambda/100/opt_acc_rate));
        mult_lambda = mult_mult * mult_lambda;
        fprintf('lambda proposal variance set to %g of previous\n\n',...
            mult_mult);
        accepted_lambda = 0 ;
        Sigma_lambda = mult_lambda * cov(log(lambda_rec(1:ii,:)));
        
        msg = fprintf('Completed: &g/%g',ii,m);
        
    end
    
    %% Update plots
    if mod(ii,upd) == 0 && doplot % Update plots every upd loops of MCMC
    
        % After the burn_in is over, we exclude it from the plots
        if ii> burn_in, startplot=burn_in; end 
        
        % If values sampled are not univariate, then cycle through cols for
        % display:
        col_theta1 = 1 + mod(col_theta1,size(theta1_rec,2));
        col_theta2 = 1 + mod(col_theta2,size(theta2_rec,2));
        col_rho    = 1 + mod(col_rho,size(rho_rec,2));
        col_lambda = 1 + mod(col_lambda,size(lambda_rec,2));
        
        % Plot
        subplot(2,2,1);
        plot(theta1_rec(startplot:ii,col_theta1),'ko');
        subplot(2,2,2);
        plot(theta2_rec(startplot:ii,col_theta2),'ko');
        subplot(2,2,3);
        plot(rho_rec(startplot:ii,col_rho),'ko');
        subplot(2,2,4);
        plot(lambda_rec(startplot:ii,col_lambda),'ko');
        
    end
    
    %% Prepare for next loop
    fprintf(repmat('/b',1,msg));
    msg = fprintf('Completed: &g/%g',ii,m);
    
end

%% Pack up and leave
theta1_os = theta1_rec .* range_t1 + min_t1;
theta2_os = theta2_rec .* range_t2 + min_t2;
results = struct('theta1',theta1_os,...
    'theta2',theta2_os,...
    'rho',rho_rec,...
    'lambda',lambda_rec,...
    'settings',settings);

end