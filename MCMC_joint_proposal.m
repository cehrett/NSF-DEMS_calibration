function [samples,Sigma] = MCMC_joint_proposal(M,burn_in,sim_xt,eta,...
    obs_x,y,Sigma_y,out_of_range,init_theta,omega,rho,lambda,...
    proposal,nugsize)
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
m = size(eta,1);
z = [y ; eta ] ; 
sim_x = sim_xt(:,1:num_cntrl) ;
sim_t = sim_xt(:,num_cntrl+1:end) ;
% The cov matrix we need, Sigma_z, can be decomposed so that this big part
% of it remains unchanged, so that we need calculate that only this once.
% Massive computation savings over getting Sigma_z from scratch each time:
Sigma_eta_xx = gp_cov(omega,sim_x,sim_x,rho,sim_t,sim_t,lambda,false);
prop_density = proposal.density ; Sigma = proposal.Sigma;

%% Initialize some variables for later use
out_of_range_rec = 0 ;
reject_rec = 0;
startplot = 1;
accepted = 0;
theta=init_theta;
msg=0;

%% Get initial log likelihood
% Set new observation input matrix:
obs_theta = repmat(theta,n,1) ; 
% Get Sigma_eta_yy:
Sigma_eta_yy = gp_cov(omega,obs_x,obs_x,rho,...
    obs_theta,obs_theta,lambda,false);
% Get Sigma_eta_xy, and hence Sigma_eta_yx
Sigma_eta_xy = gp_cov(omega,sim_x,obs_x,...
    rho,sim_t,obs_theta,lambda,false);
Sigma_eta_yx = Sigma_eta_xy';
% Combine these to get Sigma_z
Sigma_z = [ Sigma_eta_yy + Sigma_y  Sigma_eta_yx ; Sigma_eta_xy ...
    Sigma_eta_xx ] ;
% Add a nugget to Sigma_z for computational stability
Sigma_z = Sigma_z + eye(size(Sigma_z)) * nugsize(Sigma_z);
% Get log likelihood of theta_s
loglik_theta = logmvnpdf(z',0,Sigma_z) ;


%% Begin MCMC routine
for ii = 1:M
    
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
        Sigma_eta_yy = gp_cov(omega,obs_x,obs_x,rho,...
            obs_theta,obs_theta,lambda,false);
        % Get new Sigma_eta_xy, and hence Sigma_eta_yx
        Sigma_eta_xy = gp_cov(omega,sim_x,obs_x,...
            rho,sim_t,obs_theta,lambda,false);
        Sigma_eta_yx = Sigma_eta_xy';
        % Combine these to get new Sigma_z
        Sigma_z = [ Sigma_eta_yy + Sigma_y  Sigma_eta_yx ; Sigma_eta_xy ...
            Sigma_eta_xx ] ;
        % Add a nugget to Sigma_z for computational stability
        Sigma_z = Sigma_z + eye(size(Sigma_z)) * nugsize(Sigma_z);
        % Get log likelihood of theta_s
        loglik_theta_s = logmvnpdf(z',0,Sigma_z) ;
        
    end
    
    %% Get acceptance ratio statistic
    log_alpha = loglik_theta_s - loglik_theta ; 
    
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
    end
    
    %% Recordkeeping
    samples(ii,:) = theta;
    
    %% Output to console to let us know progress
    if mod(ii,10) == 0 
        fprintf(repmat('\b',1,msg));
        msg = fprintf('Completed: %g/%g\n',ii,M);
        subplot(1,2,1);
        plot(samples(startplot:end,1),'ko');
        subplot(1,2,2);
        plot(samples(startplot:end,2),'ko');
        drawnow
    end
    
    %% Tune adaptive proposal variance 
    if mod(ii,100) == 0 && ii < burn_in 
        for jj = 1 : length(theta) 
            if out_of_range_rec(jj) >= cutoff
                Sigma(jj,jj) = .75 * Sigma(jj,jj) ;
                fprintf(repmat('\b',1,msg));
                fprintf('%99.0g proposal variance reduced to %g\n',...
                    Sigma(jj,jj));
                msg = fprintf('Completed: %g/%g\n',ii,M);
            end
        end
        if all(out_of_range_rec < cutoff) && accepted < 20
            Sigma = Sigma * 0.75;
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances reduced to %g,%g\n',diag(Sigma));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        if accepted > 25
            Sigma = Sigma * 1.25;
            fprintf(repmat('\b',1,msg));
            fprintf('Proposal variances increased to %g,%g\n',diag(Sigma));
            msg = fprintf('Completed: %g/%g\n',ii,M);
        end
        accepted        = 0;
        out_of_range_rec = 0 * out_of_range_rec;
    end
    
    %% Stop plotting burn_in
    if ii > burn_in
        startplot=burn_in;
    end
    
end



end