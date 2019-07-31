function results = MCMC_discrepancy(settings)
% M is the total number of draws (inluding burn-in)
% burn_in is the length of the burn-in; may be set either as a count or as
% a proportion of M
% sim_xt is a matrix the first columns of which are the control settings 
% for the simulation observations, and the other columns of which are the
% calibration settings for those observations.
% eta is a column vector of the output of the simulation
% obs_x is the control settings for the field observations
% y is a column vector of the field observations
% Sigma_y is the observation variance
% out_of_range is a function which takes as input a proposed theta_s
% vector, and returns a logical vector describing which elements of theta_s
% are within the support of the prior on theta_s. Note that it is
% implicitly assumed here that the prior on theta_s is uniform;
% out_of_range thus defines the range of the uniform prior on theta.
% init_theta is an initial value for theta in the MCMC routine.
% proposal is a function which takes as input theta and returns a proposal
% for a new draw theta_s.
% nugsize is a function which takes as input a square matrix A and returns
% a value n to use as a nugget to stabilize A computationally.

% Define logit function for later use
logit = @(x) log(x./(1-x));

%% Unpack struct
M                       = settings.M; 
burn_in                 = settings.burn_in; 
sim_xt                  = settings.sim_xt;
eta                     = settings.eta; 
obs_x                   = settings.obs_x; 
y                       = settings.y; 
sigma2                  = settings.sigma2; 
log_sigma2_prior        = settings.log_sigma2_prior; 
log_omega_delta_prior   = settings.log_omega_delta_prior;
log_lambda_delta_prior  = settings.log_lambda_delta_prior;
out_of_range            = settings.out_of_range; 
init_theta              = settings.init_theta;
omega                   = settings.omega; 
rho                     = settings.rho; 
lambda                  = settings.lambda;
omega_delta             = settings.omega_delta_init;
lambda_delta            = settings.lambda_delta_init;
sigma2_prop_density     = settings.proposal.sigma2_prop_density ;
Sigma_sig               = settings.proposal.Sigma_sig;%pv for sigma2
omega_prop_density      = settings.proposal.omega_prop_density;
lambda_prop_density     = settings.proposal.lambda_prop_density;
prop_density            = settings.proposal.density ; %prop for theta
Sigma                   = settings.proposal.Sigma;    %pv for theta
Sigma_od                = settings.proposal.Sigma_od; %pv for omega_delta
Sigma_ld                = settings.proposal.Sigma_ld; %pv for lambda_delta
log_mh_correction_od    = settings.proposal.log_mh_correction_od;
log_mh_correction_ld    = settings.proposal.log_mh_correction_ld;
nugsize                 = settings.nugsize; 
num_out                 = settings.num_out; 
log_sig_mh_correction   = settings.log_sig_mh_correction;
log_mh_correction       = settings.log_mh_correction;
doplot                  = settings.doplot; 
log_theta_prior         = settings.log_theta_prior;
Cost_lambda             = settings.Cost_lambda;
which_outputs           = settings.which_outputs;
%X%sigma2_divs            = settings.init_sigma2_divs;
%X%sigma2_weights         =[sigma2_divs(1) sigma2_divs(2)-sigma2_divs(1)...
%X%    1-sigma2_divs(2) ] ;


%% Set plot label values
labs = [ '\omega_1' ; '\omega_2' ; '\omega_3' ; '\lambda ' ] ; 
output_labs = { 'Deflection' ; 'Rotation' ; 'Cost' };
for ii = 1:length(which_outputs)
    if which_outputs(ii) == 0
        labs(ii,:)=[]; 
    end
end

%% If burn_in is a proportion rather than a count, then convert to count
if 0<burn_in && burn_in < 1
    burn_in = ceil(burn_in * M) ;
end

%% Set cutoff value for adaptive variance
% More proposals out of support than this value (%) will cause adaptive
% variance reduction.
cutoff = 40;

%% Get some values and transformations to use later
num_cntrl = size(obs_x,2) ;
num_calib = length(init_theta) ;
n = size(obs_x,1);
num_obs = n/num_out; % Gets number of multivariate observations.
m = size(eta,1);
z = [y ; eta ] ; 
sim_x = sim_xt(:,1:num_cntrl) ;
sim_t = sim_xt(:,num_cntrl+1:end) ;
% The cov matrix we need, Sigma_z, can be decomposed so that this big part
% of it remains unchanged, so that we need calculate that only this once.
% Massive computation savings over getting Sigma_z from scratch each time:
Sigma_eta_xx = gp_cov_old(omega,sim_x,sim_x,rho,sim_t,sim_t,lambda,false);
% Get computationally nice versions of it and its inverse for post preds:
Sxx=Sigma_eta_xx + eye(size(Sigma_eta_xx)) * nugsize(Sigma_eta_xx);
iSxx = inv(Sxx);
% Now make sigma2 into a covariance matrix:
sigma2_long = repelem(sigma2,num_obs);
Sigma_y = diag(sigma2_long);

%% Initialize some variables for later use
out_of_range_rec   = zeros(size(init_theta))      ; 
reject_rec         = 0                            ;
startplot          = 10                           ;
accepted           = 0                            ;
accepted_sig       = 0                            ; % sigma2 accepted
accepted_od        = 0                            ; % omega_delta accepted
accepted_ld        = 0                            ; % lambda_delta accepted
theta              = init_theta                   ;
msg                = 0                            ;
samples            = init_theta                   ;
delta_rec          = [omega_delta lambda_delta]   ;
out_of_range_sig   = zeros(1,num_out)             ;
mult               = 10                           ; % multiplier, prop dens
mult_od            = 10                           ; % mult for pd of om_del
mult_ld            = 10                           ; % mult for pd of lam_dl
by_sample_est      = []                           ; % Collect post. preds
sds_by_sample_est  = []                           ; % Collect post sds pred

model_output.by_sample_est = by_sample_est        ; 
model_output.sds_by_sample_est = sds_by_sample_est;

%% Get initial log likelihood
% Set new observation input matrix:
obs_theta = repmat(theta,n,1) ; 
% Get Sigma_eta_yy:
Sigma_eta_yy = gp_cov_old(omega,obs_x,obs_x,rho,...
    obs_theta,obs_theta,lambda,false);
% Get Sigma_eta_xy, and hence Sigma_eta_yx
Sigma_eta_xy = gp_cov_old(omega,sim_x,obs_x,...
    rho,sim_t,obs_theta,lambda,false);
Sigma_eta_yx = Sigma_eta_xy';
% Get discrepancy covariance:
Sigma_delta = gp_cov_old(omega_delta,obs_x,obs_x,0,0,0,lambda_delta,false);
% Combine these to get Sigma_z
Sigma_z = [ Sigma_eta_yy + Sigma_y + Sigma_delta     Sigma_eta_yx   ; ...
            Sigma_eta_xy                             Sigma_eta_xx ] ;
% Add a nugget to Sigma_z for computational stability
Sigma_z = Sigma_z + eye(size(Sigma_z)) * nugsize(Sigma_z);
Sigma_z_rcond_rec = rcond(Sigma_z);
% Get log likelihood of theta
L_D = logmvnpdf(z',0,Sigma_z);
loglik_theta = L_D + log_theta_prior(theta,Cost_lambda);
% Get log likelihood of sigma2
loglik_sigma2 = L_D + log_sigma2_prior(sigma2);
% Get log likelihood of omega_delta
loglik_od = L_D + log_omega_delta_prior(omega_delta);
% Get log likelihood of lambda_delta
loglik_ld = L_D + log_lambda_delta_prior(lambda_delta);

figure('rend','painters','pos',[-1.6667    3.6667  660.0000  637.3333]);
%% Begin MCMC routine
for ii = 1:M
    
    %% Draw theta
    theta_s = prop_density(theta,Sigma) ; % Get new proposal draw
    
    if any(out_of_range(theta_s))
        
        loglik_theta_s = -Inf ; % Reject this proposed draw.
        out_of_range_rec = out_of_range_rec + out_of_range(theta_s) ;
        % Keep track of how many go out of range, for adaptive variance
        reject_rec = reject_rec  + 1; % Keep track of rejections
        
    else
        % Set new observation input matrix:
        obs_theta = repmat(theta_s,n,1) ; 
        
        %% Get new Sigma_z = Sigma_eta + [Sigma_y 0 ; 0 0], in pieces
        % Get new Sigma_eta_yy:
        Sigma_eta_yy_s = gp_cov_old(omega,obs_x,obs_x,rho,...
            obs_theta,obs_theta,lambda,false);
        % Get new Sigma_eta_xy, and hence Sigma_eta_yx
        Sigma_eta_xy_s = gp_cov_old(omega,sim_x,obs_x,...
            rho,sim_t,obs_theta,lambda,false);
        Sigma_eta_yx_s = Sigma_eta_xy_s';
        % Combine these to get new Sigma_z
        Sigma_z_s=[Sigma_eta_yy_s+Sigma_y+Sigma_delta Sigma_eta_yx_s ; ...
                   Sigma_eta_xy_s                     Sigma_eta_xx ] ;
        % Add a nugget to Sigma_z for computational stability
        Sigma_z_s = Sigma_z_s + eye(size(Sigma_z_s)) * nugsize(Sigma_z_s);
        % Get log likelihood of theta_s
        L_D_s = logmvnpdf(z',0,Sigma_z_s) ;
        loglik_theta_s = L_D_s + log_theta_prior(theta_s,Cost_lambda);
        L_D = logmvnpdf(z',0,Sigma_z) ;
        loglik_theta   = L_D + log_theta_prior(theta,Cost_lambda);
    end
    
    %% Get acceptance ratio statistic
    log_alpha = loglik_theta_s - loglik_theta + ...
        log_mh_correction(theta_s,theta); 
    
    %% Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted = accepted + 1;
    else 
        accept = 0;
    end
    if accept % Set up for next time
        loglik_theta = loglik_theta_s ;
        theta = theta_s ;
        Sigma_eta_yy = Sigma_eta_yy_s;
        Sigma_eta_yx = Sigma_eta_yx_s;
        Sigma_eta_xy = Sigma_eta_xy_s;
        Sigma_z      = Sigma_z_s;
        L_D          = L_D_s;
    end
    

    %% Draw omega_delta
    omega_delta_s = omega_prop_density(omega_delta,Sigma_od) ; 
        
    

    %% Get new Sigma_z = Sigma_eta + [Sigma_y_Sigma_delta 0 ; 0 0]
    % Set new discrep covariance matrix:
    Sigma_delta_s=...
        gp_cov_old(omega_delta_s,obs_x,obs_x,0,0,0,lambda_delta,false);
    % Combine these to get new Sigma_z
    Sigma_z_s=[Sigma_eta_yy+Sigma_y+Sigma_delta_s Sigma_eta_yx ; ...
               Sigma_eta_xy                       Sigma_eta_xx ] ;
    % Add a nugget to Sigma_z for computational stability
    Sigma_z_s = Sigma_z_s + eye(size(Sigma_z_s)) * nugsize(Sigma_z_s);
    % Get log likelihood of omega_delta_s
    L_D_s = logmvnpdf(z',0,Sigma_z_s) ;
    loglik_od_s = L_D_s + log_omega_delta_prior(omega_delta_s);
    %X%L_D = logmvnpdf(z',0,Sigma_z) ;
    loglik_od   = L_D + log_omega_delta_prior(omega_delta);
    
    %% Get acceptance ratio statistic
    log_alpha = loglik_od_s - loglik_od + ...
        log_mh_correction(omega_delta_s,omega_delta);
    
    %% Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted_od = accepted_od + 1;
    else 
        accept = 0;
    end
    if accept % Set up for next time
        loglik_od   = loglik_od_s ;
        omega_delta = omega_delta_s ;
        Sigma_delta = Sigma_delta_s ;
        Sigma_z     = Sigma_z_s;
        L_D         = L_D_s; 
    end
    
    %% Draw lambda_delta
    lambda_delta_s =lambda_prop_density(lambda_delta,Sigma_ld) ; 
        
   
    %% Get new Sigma_z = Sigma_eta + [Sigma_y_Sigma_delta 0 ; 0 0]
    % Set new discrep covariance matrix:
    Sigma_delta_s=...
        gp_cov_old(omega_delta,obs_x,obs_x,0,0,0,lambda_delta_s,false);
    % Combine these to get new Sigma_z
    Sigma_z_s=[Sigma_eta_yy+Sigma_y+Sigma_delta_s Sigma_eta_yx ; ...
               Sigma_eta_xy                       Sigma_eta_xx ] ;
    % Add a nugget to Sigma_z for computational stability
    Sigma_z_s = Sigma_z_s + eye(size(Sigma_z_s)) * nugsize(Sigma_z_s);
    % Get log likelihood of omega_delta_s
    L_D_s = logmvnpdf(z',0,Sigma_z_s) ;
    loglik_ld_s = L_D_s + log_lambda_delta_prior(lambda_delta_s);
    %X%L_D = logmvnpdf(z',0,Sigma_z) ;
    loglik_ld   = L_D + log_lambda_delta_prior(lambda_delta);
    
    %% Get acceptance ratio statistic
    log_alpha = loglik_ld_s - loglik_ld + ...
        log_mh_correction_ld(lambda_delta_s,lambda_delta);
    
    %% Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept = 1;
        accepted_ld = accepted_ld + 1;
    else 
        accept = 0;
    end
    if accept % Set up for next time
        loglik_ld    = loglik_ld_s ;
        lambda_delta = lambda_delta_s ;
        Sigma_delta  = Sigma_delta_s ;
        Sigma_z      = Sigma_z_s;
        L_D          = L_D_s; 
    end
    
    %% Recordkeeping
    samples(ii+1,:) = theta;
%     omega_delta_rec(ii+1,:) = omega_delta;
%     lambda_delta_rec(ii+1,:) = lambda_delta;
    delta_rec(ii+1,:) = [omega_delta lambda_delta];
    Sigma_z_rcond_rec(ii+1) = rcond(Sigma_z);
    
    
    %% Tune adaptive proposal variance 
    if mod(ii,100) == 0 && ii <= burn_in 
        %% Monitor likelihoods 
        fprintf(repmat('\b',1,msg));
        fprintf('loglik_theta = %f\n',loglik_theta);
        fprintf('loglik_theta_s = %f\n',loglik_theta_s);
        fprintf('log_mh_correction = %f\n',...
            log_mh_correction(theta_s,theta));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        %% Tune theta proposal variance
        if accepted < 24 
            mult = max(mult*.5,mult*accepted/24);
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances decreased\n');
            fprintf('mult = %f\n',mult);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted > 24
            %Sigma = Sigma * mult;%1.25;
%             mult = 1.25 * mult
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances increased\n');
            mult = min(mult*2,mult*accepted/24); 
            fprintf('mult = %f\n',mult);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        Sigma = cov(logit(samples)) * mult 
        msg = fprintf('Completed: %g/%g\n',ii,M);
        
        %% Tune discrepancy proposal variance
        % First omega_delta
        if accepted_od < 24 
            mult_od=max(mult_od*.5,mult_od*accepted_od/24);
            fprintf(repmat('\b',1,msg));
            fprintf('Omega delta proposal variances decreased\n');
            fprintf('mult_od = %f\n',mult_od);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted_od > 24
            fprintf(repmat('\b',1,msg));
            fprintf('Omega delta proposal variances increased\n');
            mult_od = min(mult_od*2,mult_od*accepted_od/24);
            fprintf('mult_od = %f\n',mult_od);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        Sigma_od = cov(logit(delta_rec(:,1:(end-1)))) * mult_od
        msg = fprintf('Completed: %g/%g\n',ii,M);
        % Now lambda_delta
        if accepted_ld < 44 
            mult_ld =max(mult_ld*.5,mult_ld*accepted_ld/44);
            fprintf(repmat('\b',1,msg));
            fprintf('lambda delta proposal variance increased\n');
            fprintf('mult_ld = %f\n',mult_ld);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted_ld > 44
            fprintf(repmat('\b',1,msg));
            fprintf('lambda delta proposal variance increased\n');
            mult_ld = min(mult_ld*2,mult_ld*accepted_ld/44);
            fprintf('mult_ld = %f\n',mult_ld);
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        Sigma_ld = var(log(delta_rec(:,end))) * mult_ld;
        msg = fprintf('Completed: %g/%g\n',ii,M);
        
    end
    
    %% Output to console and get posterior predictions
    if mod(ii,100) == 0
        %%% Get posterior predictions
        samps = samples(ii-99:ii,:);
        emout = ...
            em_out_many(samps,settings,0,1,Sxx,iSxx,false);
        model_output.by_sample_est = ...
            [model_output.by_sample_est ; emout.output_means];
        model_output.sds_by_sample_est = ...
            [model_output.sds_by_sample_est ; emout.output_sds];
        
        %% Print info and newline
        lag = min(50,size(samples,1)-startplot);
        vf_acf = acf(samples(startplot:ii+1,1),lag); 
        vf_acf = vf_acf(lag);
        thk_acf = acf(samples(startplot:ii+1,2),lag); 
        thk_acf = thk_acf(lag);
        fprintf(repmat('\b',1,msg));
        fprintf('accepted     = %g\n',accepted)
        fprintf('accepted_od = %g\n',accepted_od)
        fprintf('accepted_ld = %g\n',accepted_ld)
        fprintf('VF acf 50    = %g\n',vf_acf)
        fprintf('Thk acf 50   = %g\n',thk_acf)
        fprintf('\n')
        msg = fprintf('Completed: %g/%g\n',ii,M);
        
        %% Reset counters
        accepted        = 0;
        accepted_od    = 0;
        accepted_ld    = 0;
        %X% out_of_range_rec = 0 * out_of_range_rec;
        %X% out_of_range_sig = 0 * out_of_range_sig;
    end
    
    
        %% Output to console and plot to let us know progress
    if mod(ii,10) == 0 && doplot == true
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        subplot(3,4,1);
        plot(samples(startplot:end,1),'ko');
        title('Volume fraction');
        subplot(3,4,2);
        plot(samples(startplot:end,2),'ko');
        title('Thickness');
        for jj = 1:sum(size(delta_rec,2))
            subplot(3,4,jj+2);
            plot(delta_rec(startplot:end,jj),'ko');
            title(['\delta: ' labs(jj,:)]);
        end
        subplot(3,4,(sum(size(delta_rec,2))+3):8);
        plot(logit(samples(startplot:end,1)),...
            logit(samples(startplot:end,2)),'ko');
        hold on
        rr = mvnrnd(mean(logit(samples(startplot:end,:))),Sigma,250);
        plot(rr(:,1),rr(:,2),'r.');
        hold off
        for jj = 1:size(model_output.by_sample_est,2)
            subplot(3,4,8+jj)
            plot(model_output.by_sample_est(startplot:end,jj),'ko');
            title(output_labs(jj,:));
        end
        drawnow
    end
    
    %% Stop plotting burn_in
    if ii > burn_in
        startplot=burn_in;
    end
    
end

%% Pack up and leave
samples_os = samples .* settings.input_calib_ranges + ...
    settings.input_calib_mins;
results = struct('samples',samples,...
    'samples_os',samples_os,...
    'delta_samps',delta_rec,...
    'Sigma',Sigma,...
    'desired_obs',settings.desired_obs,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_delta_params',mean(delta_rec(settings.burn_in:end,:)),...
    'Sigma_z_rcond',Sigma_z_rcond_rec,...
    'model_output',model_output,...
    'settings',settings);

end