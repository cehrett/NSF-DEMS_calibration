function results = MCMC_state_aware_discrep_true_fn ( settings ) 
% This function performs state-aware calibration while including a
% discrepancy function. The input settings is a structure which is itself
% output from another function.

%% Unpack settings
M = settings.number_of_iterations; % Total # of MCMC iteratns, incl burn-in
burn_in = settings.burn_in; % Total # of iterations to treat as burn in
lambda_delta = settings.lambda_delta ; % Marginal precision of discrepancy
x = settings.cntrl_input ; % Control input points at which y is observed
xx = settings.cntrl_input_with_dum_vars; % x plus dum vars
desired_obs = settings.desired_obs ; % Target "observations", orig scale
sigma2 = settings.nugget; % Known variance of observations
mu_theta = settings.mu_theta ; % mean of GP for state-aware theta1 
dim_theta1 = settings.dim_theta1; % Dimension of state-aware parameter
dim_theta2 = settings.dim_theta2; % Dim of non-state-aware parameter
dim_nu_theta = settings.dim_nu_theta; % Dim of inputs
dim_output = settings.dim_output; % Dim of output
a_theta = settings.lambda_theta_hypers(1); % Hyperparam, shape for gamma pr
b_theta = settings.lambda_theta_hypers(2); % Hyperparam, rate for gamma pr
Sigma_xi = settings.Sigma_xi; % Covariance for proposal dist for xi
Sigma_nu_theta = settings.Sigma_nu_theta; % Cov for prop dist for nu_theta
Sigma_nu_delta = settings.Sigma_nu_delta; % Cov for prop dist for nu_delta
c_theta = settings.nu_theta_prior_param; % exp(-exp(nu_theta)) ~ Beta(1,c)
c_delta = settings.nu_delta_prior_param; % exp(-exp(nu_delta)) ~ Beta(1,c)
eta = settings.model; % The function describing the computer model
output_mean = settings.output_mean; % Prior mean of model outputs
output_sd = settings.output_sd; % Prior sd of model outputs
input_cntrl_min = settings.input_cntrl_min; % Minimum of control inputs
input_cntrl_range = settings.input_cntrl_range; % Range of control inputs
input_calib_min = settings.input_calib_min; % Minimum of calibration inputs
input_calib_range = settings.input_calib_range; % Range of calib inputs
link_fn = settings.link_fn; % Link function for the dist of theta1
which_sa = settings.which_sa; % Tells which calib inputs are state-aware

%% Standardize outputs  
y = (desired_obs - output_mean) ./ output_sd ; 
% y = repelem(y(:),numel(x),1);

%% Define logit transform and reverse transform
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x) ./ (1+exp(x));

%% Define proposal fn for theta1
theta1_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
theta2_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
% Since this proposal is not symmetric, we need to use full
% Metropolis-Hastings rather than just Metropolis. So here is the log MH
% correction for the lack of symmetry.
theta1_prop_log_mh_correction = ...
    @(t_s,t) sum(log(t_s)+log(1-t_s)-log(t)-log(1-t));


%% Normalize inputs
x = (x - input_cntrl_min) ./ input_cntrl_range;


%% Create (known) observation variance matrix
Sigma = sigma2 * eye(length(y));

%% Create receptacles for sample draws
theta1_draws       = nan(M,size(x,1),dim_theta1) ; 
xi_draws           = nan(M,dim_theta2)           ;
nu_theta_draws     = nan(M,dim_nu_theta)         ; 
lambda_theta_draws = nan(M,1)                    ;
nu_delta_draws     = nan(M,size(x,2))            ; 

%% Get initial draws of xi, nu_theta, lambda_theta, nu_delta
theta2_init             = rand(dim_theta2,1)               ;
xi_draws(1,:)           = log(-log(theta2_init))           ;
rho_theta_init          = betarnd(1,c_theta,dim_nu_theta,1);
nu_theta_draws(1,:)     = log(-log(rho_theta_init))        ;
lambda_theta_draws(1,:) = gamrnd(25,1/25)                  ;
rho_delta_init          = betarnd(1,c_delta,size(x,2),1)   ; 
nu_delta_draws(1,:)     = log(-log(rho_delta_init))        ;

%% Get initial covariance matrices
% Although gp_cov accepts lambda as an input, here we give it lambda=1, so
% that we can handle lambda elsewhere in the calculations. (Giving gp_cov
% lambda=1 is equivalent to removing lambda from gp_cov.)
R_delta=gp_cov(exp(-exp(nu_delta_draws(1,:))),x,x,1,false) + ...
    eye(size(x,1)) * 1e-4;
R_nu   =gp_cov(exp(-exp(nu_theta_draws(1,:))),x,x,1,false) +...
    eye(length(x)) * 1e-4;

%% Get initial draw of theta1
theta1_draws(1,:,:) = ones(length(x),1,dim_theta1);
while any(abs(theta1_draws(1,:,:) - .5) >= .5)
    theta1_draws(1,:,:) = mvnrnd(mu_theta*ones(length(x),1),...
        R_nu/lambda_theta_draws(1,:), dim_theta1) ;
end

%% Define scalar multipliers for adaptive proposals
mult_theta1       = 1 ;
mult_xi           = 1 ;
mult_nu_theta     = 1 ;
mult_nu_delta     = 1 ;
n_adj = 100   ; % We'll adjust these once every n_adj iterations in burn_in

%% Define counts of accepted draws
% Will be used to monitor the MCMC and to tune the proposal distributions
accepted_theta1       = 0 ; 
accepted_xi           = 0 ;
accepted_nu_theta     = 0 ;
accepted_nu_delta     = 0 ;

%% Prepare to monitor the output
msg = 0 ; % Used for console output
figure('pos',[0 0  975  675]);
colors = [      0         0    1.0000 ;
    1.0000         0         0 ;
         0    1.0000         0 ;
         0         0    0.1724 ;
    1.0000    0.1034    0.7241 ;
    1.0000    0.8276         0 ;
         0    0.3448         0 ;
    0.5172    0.5172    1.0000];
%colors = [ 'b' ; 'r' ; 'g' ; 'k' ; 'y' ; 'c' ; 'm' ] ;

%% Get initial log likelihoods
% These will be updated whenever a proposed draw is accepted in MCMC.
logL_eta = logmvnpdf( y',...
    eta(x,theta1_draws(1,:,:)',exp(-exp(xi_draws(1,:))),...
      input_cntrl_min,input_cntrl_range,...
      input_calib_min,input_calib_range,...
      output_mean,output_sd,which_sa)',...
    Sigma + R_delta / lambda_delta) ;
if dim_theta1 > 0 
    logL_SA = logmvnpdf( link_fn(theta1_draws(1,:,:)),...
        mu_theta,...
        R_nu / lambda_theta_draws(1));
    logL_prior_nu_theta = sum( (c_theta-1) * ...
        log(1-exp(-exp(nu_theta_draws(1,:))))) + ...
        sum(nu_theta_draws(1,:) - exp(nu_theta_draws(1,:)));
end
logL_prior_xi = sum(xi_draws(1,:)) - sum(exp(xi_draws(1,:)));
logL_prior_nu_delta = ...
    sum( (nu_delta_draws(1,:) - exp(nu_delta_draws(1,:))) + ...
    (c_delta - 1) * (1-exp(-exp(nu_delta_draws(1,:))))) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MCMC loop 
for ii = 2:M
    
    %% Draw theta1
    if dim_theta1 > 0
%         % Get proposal covariance (up to scale) for theta1. See Brown &
%         % Atamturktur 2016 for details; they themselves are here following 
%         % Neal 1998.
%         [U,D] = eig(R_nu) ;
%         z = mvnrnd(0*ones(length(x),1),eye(length(x)));
%         theta1_s = mult_theta1 * U *sqrt(D) * z(:) + theta1_draws(ii-1,:)';
        theta1_s = theta1_proposal(theta1_draws(ii-1,:)',R_nu)';

        % Now find whether to accept or reject theta1_s. Each of the
        % log-likelihood of theta1, theta1_s is found in 2 parts which are
        % summed; this is so those parts can be used in later parts of the
        % MCMC rather than being recalculated.
%         % Make prior on theta1 unif:
%         logL_SA=0;
        logL_theta1 = logL_eta + logL_SA; % log-likelihood of theta1

        if any(theta1_s < 0) || any(theta1_s>1)
            logL_theta1_s = -Inf;
        else
            logL_eta_s = logmvnpdf( y',...
                eta(x,theta1_s,exp(-exp(xi_draws(ii-1,:))),...
                  input_cntrl_min,input_cntrl_range,...
                  input_calib_min,input_calib_range,...
                  output_mean,output_sd,which_sa)',...
                Sigma + R_delta / lambda_delta);
            % XXX: Replacing the following lines with uniform prior bc why not?
            logL_SA_s = logmvnpdf( link_fn(theta1_s'),...
                mu_theta,...
                R_nu / lambda_theta_draws(ii-1));
%             logL_SA_s = 0;
            logL_theta1_s = logL_eta_s + logL_SA_s; %log-L of theta1_s
        end

        log_alpha = logL_theta1_s - logL_theta1 + ...
            theta1_prop_log_mh_correction(theta1_s,theta1_draws(ii-1,:)'); % log acceptance prob
        if log(rand) < log_alpha % if proposal is accepted
            accepted_theta1 = accepted_theta1 + 1 ;
            theta1_draws(ii,:) = theta1_s;
            logL_eta = logL_eta_s; % So that it can be used to draw xi
            logL_SA = logL_SA_s; % So that it can be used to draw nu_theta
        else
            theta1_draws(ii,:) = theta1_draws(ii-1,:);
        end
        
        %% Draw nu_theta
        % Sample from proposal density and make corresponding cov matrix
        nu_theta_s = ...
            mvnrnd(nu_theta_draws(ii-1,:), mult_nu_theta * Sigma_nu_theta);
        R_nu_s = gp_cov(exp(-exp(nu_theta_s)),x(:,2),x(:,2),1,false) +...
            eye(size(R_nu)) * 1e-5 ;

        % Now find whether to accept or reject nu_theta_s
        logL_nu_theta = logL_SA + logL_prior_nu_theta ; % log-L of nu_theta
        
        % Make prior on theta1 unif
%         logL_SA_s = 0;
        logL_SA_s = logmvnpdf( link_fn(theta1_draws(ii,:)),...
            mu_theta,...
            R_nu_s / lambda_theta_draws(ii-1));
        logL_prior_nu_theta_s = ...
            sum((c_theta-1) * log(1-exp(-exp(nu_theta_s)))) + ...
            sum(nu_theta_s - exp(nu_theta_s));
        logL_nu_theta_s = logL_SA_s+logL_prior_nu_theta_s;%log-L nu_theta_s

        log_alpha = logL_nu_theta_s - logL_nu_theta ; % log acceptance prob
        if log(rand) < log_alpha % if proposal is accepted
            accepted_nu_theta = accepted_nu_theta + 1 ;
            nu_theta_draws(ii,:) = nu_theta_s ; 
            logL_SA = logL_SA_s ; 
            logL_prior_nu_theta = logL_prior_nu_theta_s ; 
            R_nu = R_nu_s ; 
        else
            nu_theta_draws(ii,:) = nu_theta_draws(ii-1,:);
        end

        %% Draw lambda_theta
        % Sample from conditional density
        theta1_mmu = link_fn(theta1_draws(ii,:)-mu_theta);
        lambda_theta_draws(ii) = gamrnd( a_theta + 1/2, ...
            (b_theta + theta1_mmu(:)' *inv(R_nu)* theta1_mmu(:) / 2)^(-1));
        
    end % end of if statement: if dim_theta1 > 0
    
    if dim_theta2 > 0
        %% Draw xi
        % Sample from proposal density:
        xi_s = mvnrnd(xi_draws(ii-1,:),mult_xi * Sigma_xi) ; 

        % Now find whether to accept or reject xi_s
        logL_xi = logL_eta + logL_prior_xi; % log-likelihood of xi

        logL_eta_s = logmvnpdf( y',...
            eta(x,theta1_draws(ii,:),exp(-exp(xi_s)),...
              input_cntrl_min,input_cntrl_range,...
              input_calib_min,input_calib_range,...
              output_mean,output_sd,which_sa)',...
            Sigma + R_delta / lambda_delta) ;
        logL_prior_xi_s = sum(xi_s) - sum(exp(xi_s));
        logL_xi_s = logL_eta_s + logL_prior_xi_s; % log-likelihood of xi_s

        log_alpha = logL_xi_s - logL_xi ; % log acceptance probability
        if log(rand) < log_alpha % if proposal is accepted
            accepted_xi = accepted_xi + 1 ;
            xi_draws(ii,:) = xi_s ; 
            logL_eta = logL_eta_s ;
            logL_prior_xi = logL_prior_xi_s ;
        else
            xi_draws(ii,:) = xi_draws(ii-1,:);
        end
    end
    
    %% Draw nu_delta
    % Sample from proposal density
    nu_delta_s = mvnrnd(nu_delta_draws(ii-1,:), ...
        mult_nu_delta * Sigma_nu_delta);
    % DEBUG
%     if rand < 1/1000 disp('HEY'); end
%     nu_delta_s(1:2) = 5 ; % DEBUG
    % DEBUG
    R_delta_s = gp_cov(exp(-exp(nu_delta_s)),x,x,1,false) + ...
        eye(size(x,1)) * 1e-4;
    
    % Now find whether to accept or reject nu_delta_s
    logL_nu_delta = logL_eta + logL_prior_nu_delta ;
    
    logL_eta_s = logmvnpdf( y',...
        eta(x,theta1_draws(ii,:)',exp(-exp(xi_draws(ii,:))),...
          input_cntrl_min,input_cntrl_range,...
          input_calib_min,input_calib_range,...
          output_mean,output_sd,which_sa)',...
        Sigma + R_delta_s / lambda_delta) ;
    logL_prior_nu_delta_s = sum( (nu_delta_s - exp(nu_delta_s)) + ...
        (c_delta - 1) * (1-exp(-exp(nu_delta_s)))) ; 
    logL_nu_delta_s = logL_eta_s + logL_prior_nu_delta_s ; 
    
    log_alpha = logL_nu_delta_s - logL_nu_delta ; % log acceptance prob
    if log(rand) < log_alpha % if proposal is accepted
        accepted_nu_delta = accepted_nu_delta + 1 ; 
        nu_delta_draws(ii,:) = nu_delta_s ; 
        logL_eta = logL_eta_s ; 
        logL_prior_nu_delta = logL_prior_nu_delta_s ; 
        R_delta = R_delta_s ; 
    else
        nu_delta_draws(ii,:) = nu_delta_draws(ii-1,:);
    end
    
    %% Monitor the calibration
    if mod(ii,n_adj) == 0
        % Plot theta1 draws
        if dim_theta1 > 0
            subplot(2,3,1);
            hold on;
            for jj = 1 : size(theta1_draws,2) 
                plot(theta1_draws(1:ii,jj),...
                    'Color',colors(1+mod(jj-1,size(colors,1)),:));
            end
            titlstr = sprintf('\\theta_1 draws (%d/%d accepted)',...
                accepted_theta1,n_adj);
            title(titlstr);
            hold off;
        end
        
        % Plot theta2(:,1) draws
        subplot(2,3,2);
        hold on;
        for jj = 1 : dim_theta2
            plot(exp(-exp(xi_draws(1:ii,jj))),...
                'Color',colors(jj,:));
        end
        titlstr = sprintf('\\theta_2 draws (%d/%d accepted)',...
            accepted_xi,n_adj);
        title(titlstr);
        
        % Plot rho_theta(:,1) draws
        subplot(2,3,3);
        plot(exp(-exp(nu_theta_draws(1:ii,1))));
        titlstr = sprintf('\\rho_\\theta draws (%d/%d accepted)',...
            accepted_nu_theta,n_adj);
        title(titlstr);
        
        % Plot lambda_theta draws
        subplot(2,3,4);
        plot(lambda_theta_draws(1:ii));
        titlstr = sprintf('\\lambda_\\theta draws');
        title(titlstr);
        
        % Plot rho_delta(:,1) draws
        subplot(2,3,5);
        hold on;
        for jj = 1 : size(nu_delta_draws,2)
            plot(exp(-exp(nu_delta_draws(1:ii,jj))), ...
                'Color',colors(1+mod(jj-1,size(colors,1)),:));
        end
        titlstr = sprintf('\\rho_\\delta draws (%d/%d accepted)',...
            accepted_nu_delta,n_adj);
        title(titlstr);
        hold off;
        
%         % Plot acceptance rates
%         subplot(2,3,6);
%         hold on;
%         plot(
        
        drawnow;
        
        % Output progress to console
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %d/%d',ii,M);
        
    end
    
    %% Adjust adaptive proposal distribution covariances and mu_theta
    if mod(ii,n_adj) == 0
        if ii <= burn_in
            % Adjust scalar multipliers.
            % We adjust these multipliers to encourage acceptance rates ~ 
            % 0.23 for multivariate cases, and 0.45 for univariate cases.
            mult_theta1 = max( min( ...
                mult_theta1 * accepted_theta1 / n_adj / 0.23,...
                2 * mult_theta1), 1/2 * mult_theta1) ;
            mult_xi = max( min( ...
                mult_xi * accepted_xi / n_adj / 0.45,...
                2 * mult_xi), 1/2 * mult_xi) ; 
            mult_nu_theta = max( min( ...
                mult_nu_theta * accepted_nu_theta / n_adj / 0.23,...
                2 * mult_nu_theta), 1/2 * mult_nu_theta) ;
            mult_nu_delta = max( min( ...
                mult_nu_delta * accepted_nu_delta / n_adj / 0.23,...
                2 * mult_nu_delta), 1/2 * mult_nu_delta) ;

            % Adjust covariance matrices
            % The covariance matrices are informed by the previous draws.
            Sigma_xi = cov(xi_draws(1:ii,:)) ; 
            Sigma_nu_theta = cov(nu_theta_draws(1:ii,:)) ;
            Sigma_nu_delta = cov(nu_delta_draws(1:ii,:)) ;
            
            % Adjust mu_theta
            mu_theta = mean(theta1_draws(1:(ii-1),:));
        end
        
        % Reset accepted counters
        accepted_theta1   = 0 ; 
        accepted_xi       = 0 ;
        accepted_nu_theta = 0 ;
        accepted_nu_delta = 0 ;
    end
    
    
    
    
end % End of MCMC loop

% Output completion notice
fprintf(repmat('\b',1,msg));
fprintf('Completed: %d/%d\n\n',M,M);

results = struct(...
    'settings',settings,...
    'theta1',theta1_draws,...
    'xi',xi_draws,...
    'nu_theta',nu_theta_draws,...
    'lambda_theta',lambda_theta_draws,...
    'nu_delta',nu_delta_draws);

end % End of function