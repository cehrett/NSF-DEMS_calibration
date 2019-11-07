% scr

%% Gather results for SDoE DCTO w/ estimation of obs_var, des_var
clc ; clearvars -except dpath ; close all ; 

for discrep = 0:6 % For each version of the discrepancy

    % Turn on observation discrepancy in the model iff discrep ~=0
    if discrep == 0
        obs_discrep = false;
    else 
        obs_discrep = true;
    end
    
    % Number of chains
    m = 2;

    for ii = 1:m
        % Set initial and final observation set sizes
        obs_initial_size = 0;
        obs_final_size = 20;

        % Settings
        modular = true;
        informative_targets = false;
        des_discrep = false;
        obs_discrep_use_MLEs = false;
        obs_var = 0.05; % observation error
        des_x_size = 15;
        doplot=true;
        verbose=false;
        obs_var_est = false;

        % Set number of draws, burn_in for each mcmc:
        M = 8e3; b = .75 ; burn_in = M*b;

        % Define inputs mins and ranges 
        xmin = .5;
        xrange = .5;
        t1min = 1.5;
        t1range = 3;
        t2min = 0;
        t2range = 5;

        % Observation and target discrepancy covariance hyperparameters
        obs_rho_beta_params = [8,1];
        obs_lambda_gam_params = [8,4];
        des_rho_beta_params = [2,.4];
        des_lambda_gam_params = [40,40];

        % True theta1 function
        high_theta1 = 2.25 ; low_theta1 = 1.5 ; 
%         t1_fn = @(x) high_theta1 - ...
%             (high_theta1-low_theta1) * ...
%             exp(40*((x-t2min)/t2range)-20)./...
%             (1+exp(40*((x-t2min)/t2range)-20));
        t1_fn = @(x) high_theta1 - ...
            (high_theta1-low_theta1) * ...
            exp(40*((x-t2min)/t2range)-20)./...
            (1+exp(40*((x-t2min)/t2range)-20));


        % Define comp. model & truth (version w/ nonstandardized output)
        model_fn_ns = @(x,t1,t2) dual_calib_example_fn(...
            x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,...
            0,1,0,true) ; 
        true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
            x,xmin,xrange,...
            (t1_fn(t2*t2range+t2min)-t1min)./t1range,t1min,t1range,...
            t2,t2min,t2range,0,1,discrep,true);
            

        % Get mean and std using comp. model, define standardized version
        X=lhsdesign(1000,3); 
        Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
        mean_y = mean(Y) ; std_y = std(Y) ;
        model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
        true_phenomenon = ...
            @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

        % Get initial design
        X = lhsdesign(obs_initial_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless+ ...
            randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

        % Now set desired observations
        des_x = linspace(0,1,des_x_size)' * xrange + xmin;
        if informative_targets    
            % Find true minimum at des_x
            fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
                t1_fn(t2),t1min,t1range,...
                t2,t2min,t2range,0,1,discrep,false);
            theta2 = ...
                fmincon(fn,rand*t2range + t2min,...
                [],[],[],[],t2min,t2min+t2range);
                true_min = dual_calib_example_fn(des_x,xmin,xrange,...
                    t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
                    theta2,t2min,t2range,0,1,discrep,false);
                des_y = true_min - 2*sqrt(obs_var) ; 
        else
            des_y = 0 * ones(size(des_x,1),1);
            % Get est of target error variance
            des_y_std = (des_y - mean_y)/std_y;
            des_var = (min(des_y_std)-min((Y-mean_y)/std_y))^2;
        end

        % Since we are not using emulator, empty out simulator observations
        sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

        % And get the average discrepancy value to use as the mean
        int_fn =@(x,t) true_phenomenon(x,t) - ...
            model_fn(x,(t1_fn(t*t2range+t2min)-t1min)./t1range,t) ;
        avg_disc = integral2(int_fn,0,1,0,1) ; 
        fprintf('Average observation discrepancy: %g\n',avg_disc);
        obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

        % Emulator mean
        mean_sim = @(a,b,c) model_fn(a,b,c) ; 

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,...
            des_x,des_y,'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_final_size',obs_final_size,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'true_phenomenon',true_phenomenon,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est,...
            'des_var_est',true,...
            'true_obs_var',obs_var);

        % Perform dual calibration
        % We need a loop because sometimes an ill-conditioned matrix early
        % in the burn-in makes the whole calibration fail.
%         count = 0 ; err_count = 0 ; 
%         while count == err_count
%             try
                res = MCMC_dual_calib(settings);
%             catch ME
%                 warning('Warning: calibration failed. Retrying...');
%                 err_count = err_count + 1;
%             end
%             count = count + 1;
%             if count >= 10, rethrow(ME) ; end
%         end
        
        
        sdoe_results.theta1(:,ii) = res.theta1 ;
        sdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        sdoe_results.theta2(:,ii) = res.theta2 ;
        sdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        sdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        sdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        sdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        sdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        sdoe_results.des_rho(:,ii) = res.des_rho ;
        sdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        sdoe_results.des_lambda(:,ii) = res.des_lambda ;
        sdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        sdoe_results.obs_var(:,ii) = res.obs_var ;
        sdoe_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        sdoe_results.des_var(:,ii) = res.des_var ;
        sdoe_results.des_var_hat(ii) = mean(res.des_var(burn_in:end));
        sdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        sdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        sdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        sdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        sdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        sdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        sdoe_results.settings{ii} = settings ; 
        
        %%% Now do DCTO with all observations made ahead of time
        % Get observations
        X = lhsdesign(obs_final_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Make a col vector based on true theta1
        obs_t1 = t1_fn(obs_t2);

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless + ...
            randn(obs_final_size,1) * sqrt(obs_var) * std_y;

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,des_x,des_y,...
            'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'obs_final_size',obs_final_size,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est,...
            'des_var_est',true);

        % Perform dual calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        pdoe_results.theta1(:,ii) = res.theta1 ;
        pdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        pdoe_results.theta2(:,ii) = res.theta2 ;
        pdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        pdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        pdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        pdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        pdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        pdoe_results.des_rho(:,ii) = res.des_rho ;
        pdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        pdoe_results.des_lambda(:,ii) = res.des_lambda ;
        pdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        pdoe_results.obs_var(:,ii) = res.obs_var ;
        pdoe_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        pdoe_results.des_var(:,ii) = res.des_var ;
        pdoe_results.des_var_hat(ii) = mean(res.des_var(burn_in:end));
        pdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        pdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        pdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        pdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        pdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        pdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        pdoe_results.settings{ii} = settings ; 
        
        both_results.sdoe_results = sdoe_results;
        both_results.pdoe_results = pdoe_results;
        
        save('temp','both_results');
        
        % Close windows every now and again
        if mod(ii,10) == 0, close all ; end
        
        % Update location in the loop
        fprintf('\n####################################\n');
        fprintf('COMPLETED STEP %g/%g OF DISCREP %g\n',ii,m,discrep);
        fprintf('####################################\n');
        
        
    end
    
    % Get optimal theta2 and corresponding theta1 value
    % Get and plot true theta2
    fmfn =@(z) dual_calib_example_fn(...
        .75,xmin,xrange,t1_fn(z),t1min,t1range,...
        z,t2min,t2range,0,1,discrep,false);
    % Try two different start locations
    [theta2_1, fval_1] = fmincon(fmfn,t2min+t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [theta2_2, fval_2] = fmincon(fmfn,t2min+3*t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [~,idx] = min([fval_1 fval_2]);
    theta2s = [theta2_1 theta2_2] ; theta2 = theta2s(idx) ; 
    theta1 = t1_fn(theta2);
    
    % Add true parameter values, discrepancy version to results
    sdoe_results.discrepancy = discrep;
    pdoe_results.discrepancy = discrep;
    sdoe_results.true_theta1 = theta1;
    pdoe_results.true_theta1 = theta1;
    sdoe_results.true_theta2 = theta2;
    pdoe_results.true_theta2 = theta2;
    
    % Add results to full set of results
    results{discrep+1,1} = sdoe_results;
    results{discrep+1,2} = pdoe_results;
    
    % Save
    locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-04_SDOE_results_nobs' int2str(obs_initial_size) '-'...
    int2str(obs_final_size)]);
%     save(locstr,'results');
    
    % Close all those windows
    close all ;
    
    % Update location in the loop
    fprintf('####################################\n');
    fprintf('\nCOMPLETED DISCREPANCY VERSION %g\n',discrep);
    fprintf('####################################\n');
    
end