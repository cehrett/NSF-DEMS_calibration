function results = MCMC_dual_calib(settings)
% results = MCMC_dual_calib(settings)
%
% settings is a structure generated by the file MCMC_dual_calib_settings.m.

%% Unpack settings
M                             = settings.M;
burn_in                       = settings.burn_in;
sim_x                         = settings.sim_x;
sim_t1                        = settings.sim_t1;
sim_t2                        = settings.sim_t2;
sim_y                         = settings.sim_y;
obs_x                         = settings.obs_x;
obs_t2                        = settings.obs_t2;
obs_y                         = settings.obs_y;
des_x                         = settings.des_x;
des_y                         = settings.des_y;
min_t1                        = settings.min_t1;
range_t1                      = settings.range_t1;
min_t2                        = settings.min_t2;
range_t2                      = settings.range_t2;
obs_cov_mat                   = settings.obs_cov_mat;
mean_sim                      = settings.mean_sim;
emulator_rho                  = settings.emulator_rho;
emulator_lambda               = settings.emulator_lambda;
mean_obs                      = settings.mean_obs;
obs_rho                       = settings.obs_rho_init;
obs_lambda                    = settings.obs_lambda_init;
des_rho                       = settings.des_rho_init;
des_lambda                    = settings.des_lambda_init;
rho_proposal                  = settings.rho_proposal;
lambda_proposal               = settings.lambda_proposal;
rho_prop_log_mh_correction    = settings.rho_prop_log_mh_correction;
lambda_prop_log_mh_correction = settings.lambda_prop_log_mh_correction;
obs_Sigma_rho                 = settings.obs_rho_prop_cov;
obs_Sigma_lambda              = settings.obs_lambda_prop_cov;
des_Sigma_rho                 = settings.des_rho_prop_cov;
des_Sigma_lambda              = settings.des_lambda_prop_cov;
log_rho_prior_fn              = settings.log_rho_prior;
log_lambda_prior_fn           = settings.log_lambda_prior;
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
obs_discrep                   = settings.obs_discrep;
des_discrep                   = settings.des_discrep;
modular                       = settings.modular;

%% Set some useful variables
% nugsize tells us what size nugget to add to a matrix for computational
% stability. Here, it's just a constant, but is introduced as a function so
% that it can easily be upgraded to something fancier if desired.
nugsize = @(X) 1e-4;
% D is the vector of all outputs from all sources.
D = [ sim_y ; obs_y ; des_y ] ;
% R is the vector of the simulator outputs and real observations.
R = [ sim_y ; obs_y ] ;
opt_acc_rate = 0.234 ; % Acceptance rate treated as optimal in MCMC
upd = 100 ; % Tune adaptive covariances (and update plots) every upd loops
% The next four variables will be used when the values being sampled are
% not univariate, so that the plots cycle through which element of each
% vector is displayed.
col_theta1     = 0;
col_theta2     = 0;
col_obs_rho    = 0;
col_obs_lambda = 0;
col_des_rho    = 0;
col_des_lambda = 0;
logit = @(a) log(a./(1-a)) ; 


%% Set mean functions for GPs
% Currently, constant mean zero is used for discrepancy  GP. 
% However, it is implemented here in function form so that it can be 
% changed easily in the future if desired.
mean_des = @(a) zeros(size(a,1),1); 


%% Prepare that part of covariance matrix which will not change in MCMC
% The cov matrix we need, Sigma_z, can be decomposed so that this big part
% of it remains unchanged, so that we need calculate that only this once.
% Massive computation savings over getting Sigma_z from scratch each time:
Sigma_emulator_simsim = ...
    gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],[sim_x sim_t1 sim_t2],...
    emulator_lambda,false);


%% Initialize some variables for later use
startplot          = 1                            ;
accepted_theta1    = 0                            ; % accepted theta1
accepted_theta2    = 0                            ; % accepted theta2
accepted_des_rho   = 0                            ; % accepted des_rho
accepted_des_lambda= 0                            ; % accepted des_lambda
accepted_obs_rho   = 0                            ; % accepted obs_rho
accepted_obs_lambda= 0                            ; % accepted obs_lambda
msg                = 0                            ; % For console output
theta1_rec         = zeros(M,numel(theta1))       ;
theta2_rec         = zeros(M,numel(theta2))       ;
theta1_rec(1,:)    = theta1                       ;
theta2_rec(1,:)    = theta2                       ;
des_rho_rec        = zeros(M,numel(des_rho))      ;
des_rho_rec(1,:)   = des_rho                      ;
des_lambda_rec     = zeros(M,numel(des_lambda))   ;
des_lambda_rec(1,:)= des_lambda                   ;
obs_rho_rec        = zeros(M,numel(obs_rho))      ;
obs_rho_rec(1,:)   = obs_rho                      ;
obs_lambda_rec     = zeros(M,numel(obs_lambda))   ;
obs_lambda_rec(1,:)= obs_lambda                   ;
mult_theta1        = 10                           ; % mult. for proposal
mult_theta2        = 10                           ; 
mult_obs_rho       = 10                           ; 
mult_obs_lambda    = 10                           ; 
mult_des_rho       = 10                           ; 
mult_des_lambda    = 10                           ; 


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
    gp_cov(obs_rho,[obs_x obs_t2],[obs_x obs_t2],...
        obs_lambda,false);
Sigma_obs_dscr_obsdes = ...
    gp_cov(obs_rho,[obs_x obs_t2],...
        [des_x repmat(theta2,size(des_x,1))],...
        obs_lambda,false);
Sigma_obs_dscr_desdes = ...
    gp_cov(obs_rho,[des_x repmat(theta2,size(des_x,1))],...
        [des_x repmat(theta2,size(des_x,1))],...
        obs_lambda,false);
Sigma_des_dscr_desdes = ...
    gp_cov(des_rho,des_x,des_x,des_lambda,false);
% Now to combine everything into two big covariance matrices:
Sigma_R = [... % The covariance of R
    Sigma_emulator_simsim Sigma_emulator_simobs ; 
    Sigma_emulator_simobs' ...
        Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + obs_cov_mat...
    ];
Sigma_D = [... % The covariance of D
    Sigma_R ...
        [Sigma_emulator_simdes; ...
            Sigma_emulator_obsdes+Sigma_obs_dscr_obsdes];
    Sigma_emulator_simdes' ...
        Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes' ...
        Sigma_emulator_desdes+Sigma_obs_dscr_desdes+Sigma_des_dscr_desdes];
% Add nugget for computational tractability:
Sigma_D = Sigma_D + eye(size(Sigma_D)) * nugsize(Sigma_D) ;
% Now to get the log likelihoods
mu_R = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
    mean_sim(obs_x,repmat(theta1,size(obs_x,1),1),obs_t2) + ...
        mean_obs(obs_x,obs_t2)]; % Mean of R
mu_D = [mu_R ; 
    mean_sim(des_x,repmat(theta1,size(des_x,1),1),...
        repmat(theta2,size(des_x,1),1)) + ...
        mean_obs(des_x,repmat(theta2,size(des_x,1),1)) + ...
        mean_des(des_x)]; % Mean of D
log_cond_dens_D = logmvnpdf(D',mu_D',Sigma_D);
log_theta1_prior = log_theta1_prior_fn(theta1);
log_theta2_prior = log_theta2_prior_fn(theta2);
log_des_rho_prior = log_rho_prior_fn(des_rho);
log_des_lambda_prior = log_lambda_prior_fn(des_lambda);
log_obs_rho_prior = log_rho_prior_fn(obs_rho);
log_obs_lambda_prior = log_lambda_prior_fn(obs_lambda);

%%%%%%%%%%%%%%%
%% MCMC loop %%
%%%%%%%%%%%%%%%

figure();
for ii = 2:M
    
    %% Draw new theta1
    theta1_s = theta1_proposal(theta1,Sigma_theta1);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw theta1_s, using the updated mean and
    % covariance of the GPs. Only certain parts of the big covariance
    % matrix Sigma need to be updated -- those that depend on theta1.
    % We get mu_R and Sigma_R, and mu_D and Sigma_D. mu_R and Sigma_R are
    % restricted to the non-fake data -- the simulation runs and the real
    % observations. These are used for the modularized approach. Whether or
    % not the modularized approach is used, we also need mu_D and Sigma_D,
    % which are unrestricted in that they include all data, even the fake
    % observations. 
    mu_R_s = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
        mean_sim(obs_x,repmat(theta1_s,size(obs_x,1),1),obs_t2) + ...
            mean_obs(obs_x,obs_t2)]; % Mean of R
    mu_D_s = [mu_R_s;
        mean_sim(des_x,repmat(theta1_s,size(des_x,1),1),...
            repmat(theta2,size(des_x,1),1)) + ...
            mean_obs(des_x,repmat(theta2,size(des_x,1),1)) + ...
            mean_des(des_x)]; % Mean of D
    Sigma_emulator_simobs_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],emulator_lambda,...
        false);
    Sigma_emulator_simdes_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [des_x repmat(theta1_s,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_obsobs_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        emulator_lambda,false);
    Sigma_emulator_obsdes_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1_s,size(obs_x,1),1) obs_t2],...
        [des_x repmat(theta1_s,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_desdes_s = gp_cov(emulator_rho,...
        [des_x repmat(theta1_s,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        [des_x repmat(theta1_s,size(des_x,1),1) ...
            repmat(theta2,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_R_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs_s ;
        Sigma_emulator_simobs_s' ...
            Sigma_emulator_obsobs_s + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat
    ];
    Sigma_D_s = [...
        Sigma_R_s [Sigma_emulator_simdes_s ; ...
            Sigma_emulator_obsdes_s + Sigma_obs_dscr_obsdes];
        Sigma_emulator_simdes_s' ...
            Sigma_emulator_obsdes_s' + Sigma_obs_dscr_obsdes' ...
            Sigma_emulator_desdes_s + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes...
    ];
    % Add nugget for computational tractability
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    Sigma_R_s = Sigma_R_s + eye(size(Sigma_R_s)) * nugsize(Sigma_R_s) ;
    % Now we can get the log factors of the likelihood for the new draw,
    % using both R and D
    log_cond_dens_R_s = logmvnpdf(R',mu_R_s',Sigma_R_s);
    log_cond_dens_D_s = logmvnpdf(D',mu_D_s',Sigma_D_s);
    log_theta1_prior_s = log_theta1_prior_fn(theta1_s);
    if modular
        log_lik_theta1_s = log_cond_dens_R_s + log_theta1_prior_s; 
        % We also need to get the likelihood for the old draw, since under
        % a modular approach it is not updated in the other draws
        Sigma_R = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat
        ];
        Sigma_R = Sigma_R + eye(size(Sigma_R)) * nugsize(Sigma_R) ;
        log_cond_dens_R = logmvnpdf(R',mu_R',Sigma_R);
        log_lik_theta1 = log_cond_dens_R + log_theta1_prior;
    else
        log_lik_theta1_s = log_cond_dens_D_s + log_theta1_prior_s; 
        log_lik_theta1 = log_cond_dens_D + log_theta1_prior;
    end

    % Now we can get the acceptance ratio
    log_alpha = log_lik_theta1_s - log_lik_theta1 + ...
        theta1_prop_log_mh_correction(theta1_s,theta1);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        theta1 = theta1_s;
        mu_R = mu_R_s;
        mu_D = mu_D_s;
        Sigma_emulator_simobs = Sigma_emulator_simobs_s;
        Sigma_emulator_simdes = Sigma_emulator_simdes_s;
        Sigma_emulator_obsobs = Sigma_emulator_obsobs_s;
        Sigma_emulator_obsdes = Sigma_emulator_obsdes_s;
        Sigma_emulator_desdes = Sigma_emulator_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        log_theta1_prior = log_theta1_prior_s;
        accepted_theta1 = accepted_theta1 + 1;
    end
        
    
    %% Draw new theta2
    if sum(size(theta2))>0
    theta2_s = theta2_proposal(theta2,Sigma_theta2);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw theta2_s, using the updated mean and
    % covariance of the GPs. Only certain parts of the big covariance
    % matrix Sigma need to be updated -- those that depend on theta2.
    mu_D_s = [mean_sim(sim_x,sim_t1,sim_t2) ; ...
        mean_sim(obs_x,repmat(theta1,size(obs_x,1),1),obs_t2) + ...
            mean_obs(obs_x,obs_t2);
        mean_sim(des_x,repmat(theta1,size(des_x,1),1),...
            repmat(theta2_s,size(des_x,1),1)) + ...
            mean_obs(des_x,repmat(theta2_s,size(des_x,1),1)) + ...
            mean_des(des_x)]; % Mean of D
    Sigma_emulator_simdes_s = gp_cov(emulator_rho,[sim_x sim_t1 sim_t2],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2_s,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_obsdes_s = gp_cov(emulator_rho,...
        [obs_x repmat(theta1,size(obs_x,1),1) obs_t2],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2_s,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_emulator_desdes_s = gp_cov(emulator_rho,...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2_s,size(des_x,1),1)],...
        [des_x repmat(theta1,size(des_x,1),1) ...
            repmat(theta2_s,size(des_x,1),1)],...
        emulator_lambda,false);
    Sigma_obs_dscr_obsdes_s = gp_cov(obs_rho,...
        [obs_x obs_t2],[des_x repmat(theta2_s,size(des_x,1),1)],...
        obs_lambda,false);
    Sigma_obs_dscr_desdes_s = gp_cov(obs_rho,...
        [des_x repmat(theta2_s,size(des_x,1),1)],...
        [des_x repmat(theta2_s,size(des_x,1),1)],...
        obs_lambda,false);
    Sigma_D_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ...
            Sigma_emulator_simdes_s ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes_s + Sigma_obs_dscr_obsdes_s;
        Sigma_emulator_simdes_s' ...
            Sigma_emulator_obsdes_s' + Sigma_obs_dscr_obsdes_s' ...
            Sigma_emulator_desdes_s + Sigma_obs_dscr_desdes_s + ...
                Sigma_des_dscr_desdes...
    ];
    % Add nugget for computational tractability
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu_D_s',Sigma_D_s);
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
        mu_D = mu_D_s;
        Sigma_emulator_simdes = Sigma_emulator_simdes_s;
        Sigma_emulator_obsdes = Sigma_emulator_obsdes_s;
        Sigma_emulator_desdes = Sigma_emulator_desdes_s;
        Sigma_obs_dscr_obsdes = Sigma_obs_dscr_obsdes_s;
        Sigma_obs_dscr_desdes = Sigma_obs_dscr_desdes_s;
        log_theta2_prior = log_theta2_prior_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_theta2 = accepted_theta2 + 1;
    end
    end
    
    %% Draw new obs_rho (if discrepancy used)
    if obs_discrep
    obs_rho_s = rho_proposal(obs_rho,obs_Sigma_rho);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw obs_rho_s, using the updated covariance
    % of the GPs. Only part of the big covariance
    % matrix Sigma need to be updated -- the part that depends on obs_rho.
    % We get Sigma_R and Sigma_D. Sigma_R is
    % restricted to the non-fake data -- the simulation runs and the real
    % observations. This is used for the modularized approach. Whether or
    % not the modularized approach is used, we also need Sigma_D,
    % which is unrestricted in that it includes all data, even the fake
    % observations. 
    Sigma_obs_dscr_obsobs_s = ...
        gp_cov(obs_rho_s,[obs_x obs_t2],[obs_x obs_t2],...
            obs_lambda,false);
    Sigma_obs_dscr_obsdes_s = ...
        gp_cov(obs_rho_s,[obs_x obs_t2],...
            [des_x repmat(theta2,size(des_x,1))],...
            obs_lambda,false);
    Sigma_obs_dscr_desdes_s = ...
        gp_cov(obs_rho_s,[des_x repmat(theta2,size(des_x,1))],...
            [des_x repmat(theta2,size(des_x,1))],...
            obs_lambda,false);
    Sigma_R_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs_s + ...
                obs_cov_mat
    ];
    Sigma_D_s = [...
        Sigma_R_s [Sigma_emulator_simdes ; ...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes_s];
        Sigma_emulator_simdes' ...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes_s' ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes_s + ...
                Sigma_des_dscr_desdes...
    ];
    % Add nugget for computational tractability
    Sigma_R_s = Sigma_R_s + eye(size(Sigma_R_s)) * nugsize(Sigma_R_s) ;
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_R_s = logmvnpdf(R',mu_R',Sigma_R_s);
    log_cond_dens_D_s = logmvnpdf(D',mu_D',Sigma_D_s);
    log_obs_rho_prior_s = log_rho_prior_fn(obs_rho_s);
    % Get the log likelihoods for new and old draw either modularized or no
    if modular
        log_lik_obs_rho_s = log_cond_dens_R_s + log_obs_rho_prior_s; 
        % We also need to get the likelihood for the old draw, since under
        % a modular approach it is not updated in the other draws
        Sigma_R = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat
        ];
        Sigma_R = Sigma_R + eye(size(Sigma_R)) * nugsize(Sigma_R) ;
        log_cond_dens_R = logmvnpdf(R',mu_R',Sigma_R);
        log_lik_obs_rho = log_cond_dens_R + log_obs_rho_prior;
    else
        log_lik_obs_rho_s = log_cond_dens_D_s + log_obs_rho_prior_s; 
        log_lik_obs_rho = log_cond_dens_D + log_obs_rho_prior;
    end
    % Now we can get the acceptance ratio
    log_alpha = log_lik_obs_rho_s - log_lik_obs_rho + ...
        rho_prop_log_mh_correction(obs_rho_s,obs_rho);
    
    % Now accept obs_rho_s with probability min(alpha,1)
    if log(rand) < log_alpha
        obs_rho = obs_rho_s;
        Sigma_obs_dscr_obsobs = Sigma_obs_dscr_obsobs_s;
        Sigma_obs_dscr_obsdes = Sigma_obs_dscr_obsdes_s;
        Sigma_obs_dscr_desdes = Sigma_obs_dscr_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        log_obs_rho_prior = log_obs_rho_prior_s;
        accepted_obs_rho = accepted_obs_rho + 1;
    end
    end 
    
    %% Draw new obs_lambda (if discrepancy used)
    if obs_discrep
    obs_lambda_s = lambda_proposal(obs_lambda,obs_Sigma_lambda);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw lambda_s, using the updated covariance 
    % of the GPs. Only part of the big covariance
    % matrix Sigma need be updated -- the part that depends on obs_lambda.
    Sigma_obs_dscr_obsobs_s = ...
        gp_cov(obs_rho,[obs_x obs_t2],[obs_x obs_t2],...
            obs_lambda_s,false);
    Sigma_obs_dscr_obsdes_s = ...
        gp_cov(obs_rho,[obs_x obs_t2],...
            [des_x repmat(theta2,size(des_x,1))],...
            obs_lambda_s,false);
    Sigma_obs_dscr_desdes_s = ...
        gp_cov(obs_rho,[des_x repmat(theta2,size(des_x,1))],...
            [des_x repmat(theta2,size(des_x,1))],...
            obs_lambda_s,false);
    Sigma_R_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs_s + ...
                obs_cov_mat
    ];
    Sigma_D_s = [...
        Sigma_R_s [Sigma_emulator_simdes ; ...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes_s];
        Sigma_emulator_simdes' ...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes_s' ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes_s + ...
                Sigma_des_dscr_desdes...
    ];
    % Add nugget for computational tractability
    Sigma_R_s = Sigma_R_s + eye(size(Sigma_R_s)) * nugsize(Sigma_R_s) ;
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_R_s = logmvnpdf(R',mu_R',Sigma_R_s);
    log_cond_dens_D_s = logmvnpdf(D',mu_D',Sigma_D_s);
    log_obs_lambda_prior_s = log_lambda_prior_fn(obs_lambda_s);
    % Get the log likelihoods for new and old draw either modularized or no
    if modular
        log_lik_obs_lambda_s = log_cond_dens_R_s + log_obs_lambda_prior_s; 
        % We also need to get the likelihood for the old draw, since under
        % a modular approach it is not updated in the other draws
        Sigma_R = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat
        ];
        Sigma_R = Sigma_R + eye(size(Sigma_R)) * nugsize(Sigma_R) ;
        log_cond_dens_R = logmvnpdf(R',mu_R',Sigma_R);
        log_lik_obs_lambda = log_cond_dens_R + log_obs_lambda_prior;
    else
        log_lik_obs_lambda_s = log_cond_dens_D_s + log_obs_lambda_prior_s; 
        log_lik_obs_lambda = log_cond_dens_D + log_obs_lambda_prior;
    end
    % Now we can get the acceptance ratio
    log_alpha = log_lik_obs_lambda_s - log_lik_obs_lambda + ...
        lambda_prop_log_mh_correction(obs_lambda_s,obs_lambda);
    
    % Now accept obs_lambda_s with probability min(alpha,1)
    if log(rand) < log_alpha
        obs_lambda = obs_lambda_s;
        Sigma_obs_dscr_obsobs = Sigma_obs_dscr_obsobs_s;
        Sigma_obs_dscr_obsdes = Sigma_obs_dscr_obsdes_s;
        Sigma_obs_dscr_desdes = Sigma_obs_dscr_desdes_s;
        log_obs_lambda_prior = log_obs_lambda_prior_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_obs_lambda = accepted_obs_lambda + 1;
    end
    end
    
    %% Draw new des_rho
    if des_discrep % Only draw target outcome discrepancy if des_discrep T
    des_rho_s = rho_proposal(des_rho,des_Sigma_rho);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw des_rho_s, using the updated covariance
    % of the GPs. Only part of the big covariance
    % matrix Sigma need to be updated -- the part that depends on des_rho.
    Sigma_des_dscr_desdes_s=gp_cov(des_rho_s,des_x,des_x,des_lambda,false);
    Sigma_D_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ...
            Sigma_emulator_simdes ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes;
        Sigma_emulator_simdes' ...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes' ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes_s...
    ];
    % Add nugget for computational tractability
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu_D',Sigma_D_s);
    log_des_rho_prior_s = log_rho_prior_fn(des_rho_s);
    log_lik_des_rho_s = log_cond_dens_D_s + log_des_rho_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_des_rho = log_cond_dens_D + log_des_rho_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_des_rho_s - log_lik_des_rho + ...
        rho_prop_log_mh_correction(des_rho_s,des_rho);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        des_rho = des_rho_s;
        Sigma_des_dscr_desdes = Sigma_des_dscr_desdes_s;
        log_cond_dens_D = log_cond_dens_D_s;
        log_des_rho_prior = log_des_rho_prior_s;
        accepted_des_rho = accepted_des_rho + 1;
    end
    end
    
    %% Draw new des_lambda
    if des_discrep % Only draw target outcome discrepancy if des_discrep T
    des_lambda_s = lambda_proposal(des_lambda,des_Sigma_lambda);
    
    % Get acceptance ratio alpha
    % To do this, we'll first find the updated log factors of the
    % likelihood using the new draw lambda_s, using the updated covariance 
    % of the GPs. Only part of the big covariance
    % matrix Sigma need to be updated -- the part that depends on lambda.
    Sigma_des_dscr_desdes_s=gp_cov(des_rho,des_x,des_x,des_lambda_s,false);
    Sigma_D_s = [...
        Sigma_emulator_simsim Sigma_emulator_simobs ...
            Sigma_emulator_simdes ;
        Sigma_emulator_simobs' ...
            Sigma_emulator_obsobs + Sigma_obs_dscr_obsobs + ...
                obs_cov_mat...
            Sigma_emulator_obsdes + Sigma_obs_dscr_obsdes;
        Sigma_emulator_simdes' ...
            Sigma_emulator_obsdes' + Sigma_obs_dscr_obsdes' ...
            Sigma_emulator_desdes + Sigma_obs_dscr_desdes + ...
                Sigma_des_dscr_desdes_s...
    ];
    % Add nugget for computational tractability
    Sigma_D_s = Sigma_D_s + eye(size(Sigma_D_s)) * nugsize(Sigma_D_s) ;
    % Now we can get the log factors of the likelihood for the new draw
    log_cond_dens_D_s = logmvnpdf(D',mu_D',Sigma_D_s);
    log_des_lambda_prior_s = log_lambda_prior_fn(des_lambda_s);
    log_lik_des_lambda_s = log_cond_dens_D_s + log_des_lambda_prior_s; 
    % And we get the likelihood for the old draw
    log_lik_des_lambda = log_cond_dens_D + log_des_lambda_prior;
    % Now we can get the acceptance ratio
    log_alpha = log_lik_des_lambda_s - log_lik_des_lambda + ...
        lambda_prop_log_mh_correction(des_lambda_s,des_lambda);
    
    % Now accept theta1_s with probability min(alpha,1)
    if log(rand) < log_alpha
        des_lambda = des_lambda_s;
        Sigma_des_dscr_desdes = Sigma_des_dscr_desdes_s;
        log_des_lambda_prior = log_des_lambda_prior_s;
        log_cond_dens_D = log_cond_dens_D_s;
        accepted_des_lambda = accepted_des_lambda + 1;
    end
    end
    
    %% Record draws
    theta1_rec(ii,:)     = theta1;
    theta2_rec(ii,:)     = theta2;
    obs_rho_rec(ii,:)    = obs_rho;
    obs_lambda_rec(ii,:) = obs_lambda;
    des_rho_rec(ii,:)    = des_rho;
    des_lambda_rec(ii,:) = des_lambda;
    
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
        
        % Adjust proposal covariance for obs_rho (if discrepancy used)
        if obs_discrep
        mult_mult = max(.5,min(2,accepted_obs_rho/100/opt_acc_rate));
        mult_obs_rho = mult_mult * mult_obs_rho;
        fprintf('obs_rho proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_obs_rho = 0 ;
        obs_Sigma_rho = mult_obs_rho * cov(log(obs_rho_rec(1:ii,:)));
        end

        % Adjust proposal covariance for obs_lambda (if discrepancy used)
        if obs_discrep
        mult_mult = max(.5,min(2,accepted_obs_lambda/100/opt_acc_rate));
        mult_obs_lambda = mult_mult * mult_obs_lambda;
        fprintf('obs_lambda proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_obs_lambda = 0 ;
        obs_Sigma_lambda =mult_obs_lambda*cov(log(obs_lambda_rec(1:ii,:)));
        end
        
        % Adjust proposal covariance for des_rho (if discrepancy used)
        if des_discrep
        mult_mult = max(.5,min(2,accepted_des_rho/100/opt_acc_rate));
        mult_des_rho = mult_mult * mult_des_rho;
        fprintf('des_rho proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_des_rho = 0 ;
        des_Sigma_rho = mult_des_rho * cov(log(des_rho_rec(1:ii,:)));
        end

        % Adjust proposal covariance for des_lambda (if discrepancy used)
        if des_discrep
        mult_mult = max(.5,min(2,accepted_des_lambda/100/opt_acc_rate));
        mult_des_lambda = mult_mult * mult_des_lambda;
        fprintf('des_lambda proposal variance set to %g of previous\n',...
            mult_mult);
        accepted_des_lambda = 0 ;
        des_Sigma_lambda =mult_des_lambda*cov(log(des_lambda_rec(1:ii,:)));
        end
        
        fprintf('\n');
        msg = fprintf('Completed: %g/%g  ',ii,M);
        
    end
    
    %% Update plots
    if mod(ii,upd) == 0 && doplot % Update plots every upd loops of MCMC
    
        % After the burn_in is over, we exclude it from the plots
        if ii> burn_in, startplot=burn_in; end 
        
        % If values sampled are not univariate, then cycle through cols for
        % display:
        col_theta1     = 1 + mod(col_theta1,size(theta1_rec,2));
        col_theta2     = 1 + mod(col_theta2,size(theta2_rec,2));
        col_obs_rho    = 1 + mod(col_obs_rho,size(obs_rho_rec,2));
        col_obs_lambda = 1 + mod(col_obs_lambda,size(obs_lambda_rec,2));
        col_des_rho    = 1 + mod(col_des_rho,size(des_rho_rec,2));
        col_des_lambda = 1 + mod(col_des_lambda,size(des_lambda_rec,2));
        
        % Plot
        subplot(2,3,1);
        plot(theta1_rec(startplot:ii,col_theta1),'ko');
        subplot(2,3,2);
        if prod(size(theta2_rec))>0
            plot(theta2_rec(startplot:ii,col_theta2),'ko');
        end
        subplot(2,3,3);
        plot(obs_rho_rec(startplot:ii,col_obs_rho),'ko');
        subplot(2,3,4);
        plot(obs_lambda_rec(startplot:ii,col_obs_lambda),'ko');
        % Only need plot these two if des_discrep == true
        if des_discrep
            subplot(2,3,5);
            plot(des_rho_rec(startplot:ii,col_des_rho),'ko');
            subplot(2,3,6);
            plot(des_lambda_rec(startplot:ii,col_des_lambda),'ko');
        end
        
        drawnow;
    end
    
    %% Prepare for next loop
    fprintf(repmat('\b',1,msg));
    msg = fprintf('Completed: %g/%g  ',ii,M);
    
end

%% Pack up and leave
theta1_os = theta1_rec .* range_t1 + min_t1;
if prod(size(theta2_rec))>0, theta2_os = theta2_rec .* range_t2 + min_t2;
else theta2_os = theta2_rec ; end
results = struct('theta1',theta1_os,...
    'theta2',theta2_os,...
    'obs_rho',obs_rho_rec,...
    'obs_lambda',obs_lambda_rec,...
    'des_rho',des_rho_rec,...
    'des_lambda',des_lambda_rec,...
    'settings',settings);

end