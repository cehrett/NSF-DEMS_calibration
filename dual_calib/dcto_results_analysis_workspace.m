%% Check out results from SDoE vs PDoE
clc ; clearvars -except dpath ; close all ;

obs_initial_size = 0 ; obs_final_size = 20;
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-10-01_SDOE_results']);
locstr = [dpath,'dual_calib\dual_calib_stored_data\',...
    '2019-10-04_SDOE_results_nobs' int2str(obs_initial_size) '-'...
    int2str(obs_final_size)];
load(locstr,'results');

for ii=1:7
    true_theta1 = results{ii,1}.true_theta1 ;
    sdoe_t1 = mean(results{ii,1}.theta1_hat) ;
    sdoe_t1var = var(results{ii,1}.theta1_hat) ;
    pdoe_t1 = mean(results{ii,2}.theta1_hat);
    pdoe_t1var = var(results{ii,2}.theta1_hat) ;
    fprintf(['For discrep %g, true theta1 = %g.\n'...
        'SDoE got point estimate %g with variance %g.\n'...
        'PDoE got point estimate %g with variance %g.\n\n'],...
        ii-1,true_theta1,sdoe_t1,sdoe_t1var,pdoe_t1,pdoe_t1var);
    
    true_theta2 = results{ii,1}.true_theta2 ;
    sdoe_t2 = mean(results{ii,1}.theta2_hat) ;
    sdoe_t2var = var(results{ii,1}.theta2_hat) ;
    pdoe_t2 = mean(results{ii,2}.theta2_hat);
    pdoe_t2var = var(results{ii,2}.theta2_hat) ;
    fprintf(['True theta2 = %g.\n'...
        'SDoE got point estimate %g with variance %g.\n'...
        'PDoE got point estimate %g with variance %g.\n'],...
        true_theta2,sdoe_t2,sdoe_t2var,pdoe_t2,pdoe_t2var);
    
    fprintf('######################################################\n');
    
end

%% Check out results from DCTO vs KOH+CTO
clc ; clearvars -except dpath ; close all ;

locstr = [dpath,'dual_calib\dual_calib_stored_data\',...
    '2019-10-04_DCTO_vs_KOHCTO_results'];
load(locstr);

for ii=1:7
    true_theta1 = results{ii,1}.true_theta1 ;
    dcto_t1 = mean(results{ii,1}.theta1_hat) ;
    dcto_t1var = var(results{ii,1}.theta1_hat) ;
    khtc_t1 = mean(results{ii,2}.theta1_hat);
    khtc_t1var = var(results{ii,2}.theta1_hat) ;
    fprintf(['For discrep %g, true theta1 = %g.\n'...
        'DCTO got point estimate %g with variance %g.\n'...
        'KOHC got point estimate %g with variance %g.\n\n'],...
        ii-1,true_theta1,dcto_t1,dcto_t1var,khtc_t1,khtc_t1var);
    
    true_theta2 = results{ii,1}.true_theta2 ;
    dcto_t2 = mean(results{ii,1}.theta2_hat) ;
    dcto_t2var = var(results{ii,1}.theta2_hat) ;
    ctod_t2 = mean(results{ii,3}.theta2_hat);
    ctod_t2var = var(results{ii,3}.theta2_hat) ;
    fprintf(['True theta2 = %g.\n'...
        'DCTO got point estimate %g with variance %g.\n'...
        'CTOD got point estimate %g with variance %g.\n'],...
        true_theta2,dcto_t2,dcto_t2var,ctod_t2,ctod_t2var);
    
    fprintf('######################################################\n');
    
end

%% Check out results from DCTO vs KOH+CTO w/ estimated des_var
clc ; clearvars -except dpath ; close all ;

% locstr = [dpath,...
%     'dual_calib\dual_calib_stored_data\',...
%     '2019-11-20_DCTO_vs_KOHCTO_results'];
locstr = [dpath,...
    'dual_calib\dual_calib_stored_data\',...
    '2020-07-18_DCTO_vs_KOHCTO_results'];
load(locstr);

% Helper function
rmse = @(truval,samps) sqrt( sum((samps-truval).^2)/length(samps) ) ;

for ii=1:7
    true_theta1 = results{ii,1}.true_theta1 ;
    dcto_t1 = mean(results{ii,1}.theta1_hat) ;
%     dcto_t1var = var(results{ii,1}.theta1_hat) ;
    dcto_t1var = mean(results{ii,1}.theta1_var);
    dcto_t1rmse = rmse(true_theta1,results{ii,1}.theta1_hat);
    khtc_t1 = mean(results{ii,2}.theta1_hat);
%     khtc_t1var = var(results{ii,2}.theta1_hat) ;
    khtc_t1var = mean(results{ii,2}.theta1_var);
    khtc_t1rmse = rmse(true_theta1,results{ii,2}.theta1_hat);
    fprintf(['For discrep %g, true theta1 = %3.3g.\n'...
        'DCTO got point estimate %3.3g with RMSE %3.3g\n'...
        '                              and var. %3.3g.\n'...
        'KOHC got point estimate %3.3g with RMSE %3.3g\n'...
        '                              and var. %3.3g.\n\n'],...
        ii-1,true_theta1,dcto_t1,dcto_t1rmse,dcto_t1var,...
        khtc_t1,khtc_t1rmse,khtc_t1var);
    
    true_theta2 = results{ii,1}.true_theta2 ;
    dcto_t2 = mean(results{ii,1}.theta2_hat) ;
%     dcto_t2var = var(results{ii,1}.theta2_hat) ;
    dcto_t2var = mean(results{ii,1}.theta2_var);
    dcto_t2rmse = rmse(true_theta2,results{ii,1}.theta2_hat);
    ctod_t2 = mean(results{ii,3}.theta2_hat);
%     ctod_t2var = var(results{ii,3}.theta2_hat) ;
    ctod_t2var = mean(results{ii,3}.theta2_var);
    ctod_t2rmse = rmse(true_theta2,results{ii,3}.theta2_hat);
    fprintf(['True theta2 = %3.3g.\n'...
        'DCTO got point estimate %3.3g with RMSE %3.3g\n'...
        '                              and var. %3.3g.\n'...
        'CTOD got point estimate %3.3g with RMSE %3.3g\n'...
        '                              and var. %3.3g.\n\n'],...
        true_theta2,dcto_t2,dcto_t2rmse,dcto_t2var,...
        ctod_t2,ctod_t2rmse,ctod_t2var);
    
%     true_obs_var = results{ii,1}.settings{1}.obs_var ;
%     dcto_obs_var = mean(results{ii,1}.obs_var_hat) ;
%     dcto_obs_var_var = var(results{ii,1}.obs_var_hat) ;
%     khtc_obs_var = mean(results{ii,2}.obs_var_hat);
%     khtc_obs_var_var = var(results{ii,2}.obs_var_hat) ;
%     fprintf(['True obs_var = %g.\n'...
%         'DCTO got point estimate %g with variance %g.\n'...
%         'CTOD got point estimate %g with variance %g.\n\n'],...
%         true_obs_var,dcto_obs_var,dcto_obs_var_var,khtc_obs_var,...
%         khtc_obs_var_var);
    
%     true_des_var = results{ii,1}.settings{1}.des_var ;
    dcto_des_var = mean(results{ii,1}.des_var_hat) ;
    dcto_des_var_var = var(results{ii,1}.des_var_hat) ;
    ctod_des_var = mean(results{ii,3}.des_var_hat);
    ctod_des_var_var = var(results{ii,3}.des_var_hat) ;
    fprintf(['True des_var = NA.\n'...
        'DCTO got point estimate %3.3g with variance %3.3g.\n'...
        'CTOD got point estimate %3.3g with variance %3.3g.\n'],...
        dcto_des_var,dcto_des_var_var,ctod_des_var,...
        ctod_des_var_var);
    
    fprintf('######################################################\n');
    
end



%% Check out results from SDoE vs PDoE w/ estimated des_var
clc ; clearvars -except dpath ; close all ;

obs_initial_size = 0 ; obs_final_size = 20;
locstr = [dpath,'dual_calib\dual_calib_stored_data\',...
    '2019-10-31_SDOE_results_desvarest_nobs' ...
    int2str(obs_initial_size) '-'...
    int2str(obs_final_size)];
load(locstr,'results');

% Helper function
rmse = @(truval,samps) sqrt( sum((samps-truval).^2)/length(samps) ) ;

for ii=1:7
    true_theta1 = results{ii,1}.true_theta1 ;
    sdoe_t1 = mean(results{ii,1}.theta1_hat) ;
    sdoe_t1var = mean(results{ii,1}.theta1_var) ;
    sdoe_t1rmse = rmse(true_theta1,results{ii,1}.theta1_hat);
    pdoe_t1 = mean(results{ii,2}.theta1_hat);
    pdoe_t1var = mean(results{ii,2}.theta1_var) ;
    pdoe_t1rmse = rmse(true_theta1,results{ii,2}.theta1_hat);
    fprintf(['For discrep %g, true theta1 = %g.\n'...
        'SDoE got point estimate %3.3g with rmse %3.3g\n'...
        '                        and variance %3.3g.\n'...
        'PDoE got point estimate %3.3g with rmse %3.3g\n'...
        '                        and variance %3.3g.\n\n'],...
        ii-1,true_theta1,sdoe_t1,sdoe_t1rmse,sdoe_t1var,...
        pdoe_t1,pdoe_t1rmse,pdoe_t1var);
    
    true_theta2 = results{ii,1}.true_theta2 ;
    sdoe_t2 = mean(results{ii,1}.theta2_hat) ;
    sdoe_t2var = var(results{ii,1}.theta2_hat) ;
    sdoe_t2rmse = rmse(true_theta2,results{ii,1}.theta2_hat);
    pdoe_t2 = mean(results{ii,2}.theta2_hat);
    pdoe_t2var = var(results{ii,2}.theta2_hat) ;
    pdoe_t2rmse = rmse(true_theta2,results{ii,2}.theta2_hat);
    fprintf(['True theta2 = %g.\n'...
        'SDoE got point estimate %3.3g with rmse %3.3g\n'...
        '                         and variance %3.3g.\n'...
        'PDoE got point estimate %3.3g with rmse %3.3g\n'...
        '                         and variance %3.3g.\n\n'],...
        true_theta2,sdoe_t2,sdoe_t2rmse,sdoe_t2var,...
        pdoe_t2,pdoe_t2rmse,pdoe_t2var);
    
    dcto_des_var = mean(results{ii,1}.des_var_hat) ;
    dcto_des_var_var = var(results{ii,1}.des_var_hat) ;
    pdoe_des_var = mean(results{ii,2}.des_var_hat);
    pdoe_des_var_var = var(results{ii,2}.des_var_hat) ;
    fprintf(['True des_var = NA.\n'...
        'DCTO got point estimate %3.3g with variance %g.\n'...
        'CTOD got point estimate %3.3g with variance %g.\n'],...
        dcto_des_var,dcto_des_var_var,pdoe_des_var,...
        pdoe_des_var_var);
    
    fprintf('######################################################\n');
    
end

%% Helper functions

