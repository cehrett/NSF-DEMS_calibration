function [objectives, cons] = emulator_mean(res, x, theta)
    % Gets emulator output at a given theta value.
    % res is a set of results from CTO. This is used for convenience, since
    % res.settings has the emulator hyperparameters we need.
    % x is a set of operational domain inputs, on original scale.
    % theta is a set of design variables at which we want the emulator
    % mean output. It should be provided on the original scale.
    % This function only outputs the mean, no
    % variance/covariance/uncertainty estimation.
    % The output is given on the original scale.
    % cons is just to make it work with ngpm package for NSGA-II.
    %
    % It should be easy to mod this to output the variance of each output
    % too.
    % This version last updated 2020-09-16.
    
    cons=[];
    
    
    % Define inputs mins and ranges
    xmin = res.settings.min_x;
    xrange = res.settings.range_x;
    sizex = size(xmin,2);
    tmin = [res.settings.min_t1 res.settings.min_t2];
    trange = [res.settings.range_t1 res.settings.range_t2];
    
    % Number of outputs
    num_out = size(res.settings.emulator_rho,2);

    
    % Get set up for loop that will find posterior predictive dist
    pred_xt_all = ([x theta] - [xmin tmin])./[xrange trange];
    sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
    prior_mean_xt_all = res.settings.mean_sim(pred_xt_all) ; 
    eta = res.settings.sim_y;
    
    prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 
    cov_xt_xt = nan(size(sim_xt,1),size(sim_xt,1),num_out);
    for ii=1:num_out
        cov_xt_xt(:,:,ii) = ...
            gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
            res.settings.emulator_lambda(ii),false);
        cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
            1e-5*eye(size(cov_xt_xt(:,:,ii)));
    end

    % Get posterior predictive distribution for each output
    % Only do maxsimul samples at once, to ease computational expense
    maxsimul = 500;
    for jj = 1 : ceil(size(theta,1)/maxsimul)
        startpt = (jj-1) * maxsimul + 1 ;
        endpt = min(jj * maxsimul, size(theta,1));
        pred_xt = pred_xt_all(startpt:endpt,:);
        for ii=1:num_out
            cov_pxt_xt = ...
                gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
                res.settings.emulator_lambda(ii),false);
            cov_pxt_pxt = ...
                gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
                res.settings.emulator_lambda(ii),false);
            pred_mean_xt_current(:,ii) = ...
                prior_mean_xt_all(startpt:endpt,ii) + ...
                cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
                (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
            pred_cov_xt_current(:,:,ii) = ...
                cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
%             fprintf('\n%g of %g, %g of 3 Done\n',jj,...
%                 ceil(size(theta,1)/maxsimul),ii);
            pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
            pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
                pred_cov_xt_current(:,:,ii);
        end
    end
    % Put back on original scale
    mean_y = res.settings.mean_y; std_y = res.settings.std_y;
    pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
    objectives = pred_mean_xt_os';
    
end