%scr
%% Check out results from SDoE vs PDoE
clc ; clearvars -except dpath ; close all ;

obs_initial_size = 0 ; obs_final_size = 20;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-04_SDOE_results_nobs' int2str(obs_initial_size) '-'...
    int2str(obs_final_size)]);
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-10-01_SDOE_results']);
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

locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-04_DCTO_vs_KOHCTO_results']);
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