function result = ...
    MCMC_set_total_obs_var_true_fn(settings)
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
out_of_range            = settings.out_of_range; 
init_theta              = settings.init_theta;
omega                   = settings.omega; 
rho                     = settings.rho; 
lambda                  = settings.lambda;
proposal                = settings.proposal; 
nugsize                 = settings.nugsize; 
num_out                 = settings.num_out; 
log_sig_mh_correction   = settings.log_sig_mh_correction;
log_mh_correction       = settings.log_mh_correction;
doplot                  = settings.doplot; 
log_theta_prior         = settings.log_theta_prior;
Cost_lambda             = settings.Cost_lambda;
which_outputs           = settings.which_outputs;
sigma2_divs             = settings.init_sigma2_divs;
sigma2_weights          = [ sigma2_divs(1) sigma2_divs(2)-sigma2_divs(1)...
                            1-sigma2_divs(2) ] ;
y_means                 = settings.output_means';
y_sds                   = settings.output_sds';
theta_ranges            = settings.input_calib_ranges;
theta_mins              = settings.input_calib_mins;


%% Set plot label values
labs = [ 'defl' ; 'rotn' ; 'cost' ] ; 
for ii = 1:length(which_outputs)
    if which_outputs(ii) == 0
        labs(ii,:)=[]; 
    end
end

%% If burn_in is a proportion rather than a count, then convert to count
if 0<burn_in && burn_in < 1
    burn_in = ceil(burn_in * M) ;
end

%% Set proposal density for sigma2:
sigma2_prop_density = proposal.sigma2_prop_density ;
Sigma_sig = proposal.Sigma_sig ; % this is the proposal var for sigma2

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
%X% z = [y ; eta ] ; 
sim_x = sim_xt(:,1:num_cntrl) ;
sim_t = sim_xt(:,num_cntrl+1:end) ;
% The cov matrix we need, Sigma_z, can be decomposed so that this big part
% of it remains unchanged, so that we need calculate that only this once.
% Massive computation savings over getting Sigma_z from scratch each time:
%X% Sigma_eta_xx = gp_cov(omega,sim_x,sim_x,rho,sim_t,sim_t,lambda,false);
prop_density = proposal.density ; 
Sigma = proposal.Sigma;
% Now make sigma2 into a covariance matrix:
sigma2_long = repelem(sigma2*sigma2_weights/sum(sigma2_weights),num_obs);
Sigma_y = diag(sigma2_long);

%% Initialize some variables for later use
out_of_range_rec   = zeros(size(init_theta)); 
reject_rec         = 0                ;
startplot          = 10               ;
accepted           = 0                ;
accepted_sig       = 0                ;
theta              = init_theta       ;
msg                = 0                ;
samples            = init_theta       ;
sigma2_weights_rec = sigma2_weights   ;
out_of_range_sig   = zeros(1,num_out) ;
mult               = 10               ;

%% Get initial log likelihood
% Set new observation input matrix:
obs_theta = repmat(theta,n,1) ; 
% Get Sigma_eta_yy:
%X% Sigma_eta_yy = gp_cov(omega,obs_x,obs_x,rho,...
%X%    obs_theta,obs_theta,lambda,false);
% Get Sigma_eta_xy, and hence Sigma_eta_yx
%X%Sigma_eta_xy = gp_cov(omega,sim_x,obs_x,...
%X%    rho,sim_t,obs_theta,lambda,false);
%X%Sigma_eta_yx = Sigma_eta_xy';
% Combine these to get Sigma_z
%X% Sigma_z = [ Sigma_eta_yy + Sigma_y  Sigma_eta_yx ; Sigma_eta_xy ...
%X%    Sigma_eta_xx ] ;
% Add a nugget to Sigma_z for computational stability
%X% Sigma_z = Sigma_z + eye(size(Sigma_z)) * nugsize(Sigma_z);
% Get log likelihood of theta
% Note that I here assume that the desired observations are, in order, for
% some k, k observations of deflection, then k observations of rotation at
% those same settings, then k observations of cost at those same settings
% again.
obs_x_os = obs_x(1:size(obs_x,1)/3,end) ...
    * settings.input_cntrl_ranges + settings.input_cntrl_mins;
theta_os = theta .* theta_ranges + theta_mins ;
obs_x_theta = [obs_x_os repmat(theta_os,size(obs_x_os,1),1)];
true_y = (Ex_sim(obs_x_theta) - y_means)./y_sds;
true_y = true_y(:);
L_D = logmvnpdf(y',true_y,Sigma_y);
loglik_theta = L_D + log_theta_prior(theta,Cost_lambda);
% Get log likelihood of sigma2
loglik_sigma2 = L_D + log_sigma2_prior(sigma2);

if doplot figure(); end % For observing MCMC
%% Begin MCMC routine
for ii = 1:M
    
    theta_s = prop_density(theta,Sigma) ; % Get new proposal draw
    theta_s_os = theta_s .* theta_ranges + theta_mins;
    
    if any(out_of_range(theta_s))
        
        loglik_theta_s = -Inf ; % Reject this proposed draw.
        out_of_range_rec = out_of_range_rec + out_of_range(theta_s) ;
        % Keep track of how many go out of range, for adaptive variance
        reject_rec = reject_rec  + 1; % Keep track of rejections
        
    else
        
        obs_x_theta_s = [obs_x_os repmat(theta_s_os,size(obs_x_os,1),1)];
        true_y_s = (Ex_sim(obs_x_theta_s) - y_means)./y_sds;
        true_y_s = true_y_s(:);
        L_D_s = logmvnpdf(y',true_y_s',Sigma_y);
        loglik_theta_s = L_D_s + log_theta_prior(theta_s,Cost_lambda);
        L_D = logmvnpdf(y',true_y',Sigma_y) ;
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
        true_y = true_y_s;
        obs_x_theta = obs_x_theta_s;
    end
    
    %% Propose new sigma2 weights
    sigma2_divs_s = sigma2_prop_density(sigma2_divs,Sigma_sig);
    sigma2_weights_s = [ sigma2_divs_s(1) ...
        sigma2_divs_s(2)-sigma2_divs_s(1) 1-sigma2_divs_s(2) ] ;
    
    %% Get loglikelihood
    if any([sigma2_weights_s sigma2_divs_s] <= 0)
        out_of_range_sig = out_of_range_sig + (sigma2_weights_s<=0) + ...
            (sigma2_weights_s>=1);
        loglik_sigma2_s = -Inf;
    else
        %% Now make sigma2_s into a covariance matrix:
        sigma2_props_s = sigma2_weights_s/sum(sigma2_weights_s);
        %Xd sigma_2_s = sigma2_weights_s .* sigma_2;
        sigma2_s_long = repelem(sigma2 * sigma2_props_s,num_obs);
        Sigma_y_s = diag(sigma2_s_long);

        
        loglik_sigma2_s = logmvnpdf(y',true_y',Sigma_y_s) + ...
            log_sigma2_prior(sigma2_weights_s);
        loglik_sigma2   = logmvnpdf(y',true_y',Sigma_y) + ...
            log_sigma2_prior(sigma2_weights);
    end
    
    %% Get acceptance ratio statistic
    log_alpha = loglik_sigma2_s - loglik_sigma2 + ...
        log_sig_mh_correction(sigma2_weights_s,sigma2_weights) ; 
    
    %% Randomly accept or reject with prob. alpha; update accordingly
    if log(rand) < log_alpha
        accept_sig = 1;
        accepted_sig = accepted_sig + 1;
    else 
        accept_sig = 0;
    end
    if accept_sig % Set up for next time
        loglik_sigma2 = loglik_sigma2_s ;
        sigma2_weights = sigma2_weights_s ;
        sigma2_divs    = sigma2_divs_s ;
        Sigma_y = Sigma_y_s;
    end
    
    %% Recordkeeping
    samples(ii+1,:) = theta;
    sigma2_weights_rec(ii+1,:) = sigma2_weights;
    
    %% Tune adaptive proposal variance 
    if mod(ii,100) == 0 && ii <= burn_in 
        %% Tune theta proposal variance
        if accepted < 24 mult = max(mult*.75,mult*accepted/24)
        end
        if accepted > 24
            %Sigma = Sigma * mult;%1.25;
%             mult = 1.25 * mult
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances increased\n');
            mult = min(mult*2,mult*accepted/24) 
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        Sigma = cov(logit(samples)) * mult 
        msg = fprintf('Completed: %g/%g\n',ii,M);
        
        %% Tune sigma2 proposal variance
        if accepted_sig < 24
            mult_sig = max(0.75,accepted_sig/24);
            Sigma_sig = Sigma_sig * mult_sig;
            fprintf(repmat('\b',1,msg));
            fprintf('sigma2 proposal variances reduced to %g\n',...
                diag(Sigma_sig))
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted_sig > 24
            mult_sig = min(2,accepted_sig/24);
            Sigma_sig = Sigma_sig * mult_sig;
            fprintf(repmat('\b',1,msg));
            fprintf('sigma2 proposal variances increased to %g\n',...
                diag(Sigma_sig))
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        
    end
    
    if mod(ii,100) == 0
        % Print info and newline
        lag = min(50,size(samples,1)-startplot);
        vf_acf = acf(samples(startplot:ii+1,1),lag); 
        vf_acf = vf_acf(lag);
        thk_acf = acf(samples(startplot:ii+1,2),lag); 
        thk_acf = thk_acf(lag);
        fprintf(repmat('\b',1,msg));
        fprintf('accepted     = %g\n',accepted)
        fprintf('accepted_sig = %g\n',accepted_sig)
        fprintf('VF acf 50    = %g\n',vf_acf)
        fprintf('Thk acf 50   = %g\n',thk_acf)
        fprintf('\n')
        msg = fprintf('Completed: %g/%g\n',ii,M);
        
        %% Reset counters
        accepted        = 0;
        accepted_sig    = 0;
        out_of_range_rec = 0 * out_of_range_rec;
        out_of_range_sig = 0 * out_of_range_sig;
    end
    
    
        %% Output to console to let us know progress
    if mod(ii,10) == 0 && doplot == true
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        subplot(2,3,1);
        plot(samples(startplot:end,1),'ko');
        title('Volume fraction');
        subplot(2,3,2);
        plot(samples(startplot:end,2),'ko');
        title('Thickness');
        for jj = 1:num_out
            subplot(2,3,jj+2);
            % Comment out old version
%             plot(sigma2_weights_rec(startplot:end,jj)./...
%                 sum(sigma2_weights_rec(startplot:end,:),2),'ko');
            plot(sigma2_weights_rec(startplot:end,jj),'ko');
            title(['\sigma^2: ' labs(jj,:)]);
        end
        subplot(2,3,6);
        plot(logit(samples(startplot:end,1)),...
            logit(samples(startplot:end,2)),'ko');
        hold on
        rr = mvnrnd(mean(logit(samples(startplot:end,:))),Sigma,200);
        plot(rr(:,1),rr(:,2),'r.');
        hold off
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
result = struct('samples',samples,...
    'samples_os',samples_os,...
    'sigma2_weights',sigma2_weights_rec,...
    'Sigma',Sigma,...
    'desired_obs',settings.desired_obs,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_sigma2',mean(sigma2_weights_rec(settings.burn_in:end,:)),...
    'settings',settings);

end