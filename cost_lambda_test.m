clc; clear all; close all;

%% Add paths
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');


%% User defined values
M = 1e4;
desired_obs = [.65 0.77];
which_outputs = [ 1 1 0 ] ; %Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
results2 = cell(3,1);
count = 1;

for ii = [50]
    settings.Cost_lambda=ii;
    
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);


    post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
        settings.sim_xt,settings.eta,settings.output_sds,...
        settings.output_means,settings.omega,settings.rho,...
        settings.lambda,which_outputs);
    result = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
        'init',samples(1,:),'desired_data',desired_obs,...
        'sigma2_prior',settings.log_sigma2_prior,...
        'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
        'proposal',settings.proposal,'nugsize',settings.nugsize,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'Cost_lambda',settings.Cost_lambda);
    
    results2{count}=result;
    count = count + 1;
end

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
'results_Cost_lambda_test'],...
'results');

%% Figures

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
'results_Cost_lambda_test'],...
'results');
r=cell(4,1); r{1}=results{1}; r{2}=results{10};r{3}=results2{1};r{4}=results{100};results=r;
burn_in=2002; logit = @(x) log(x./(1-x));

h=figure('rend','painters','pos',[10 10 800 600]);
for ii =1:4
    
    subplot(2,3,1);
    plot(results{ii}.samples(burn_in:end,1),'ko');
    title('Volume fraction');
    subplot(2,3,2);
    plot(results{ii}.samples(burn_in:end,2),'ko');
    title('Thickness');
    subplot(2,3,4);
    plot(results{ii}.sigma2(burn_in:end,1),'ko');
    title('Deflection \sigma^2');
    subplot(2,3,5);
    plot(results{ii}.sigma2(burn_in:end,2),'ko');
    title('Rotation \sigma^2');

    subplot(2,3,6);
    lbsamps = logit(results{ii}.samples(burn_in:end,:));
    plot(lbsamps(:,1),lbsamps(:,2),'ko');
    hold on;
    rr = mvnrnd(mean(lbsamps),results{ii}.Sigma,150);
    plot(rr(:,1),rr(:,2),'r.');
    hold off;
    title('(Logit) samples with proposal cov.');

    suptitle(strcat('Prior on VF and Thickness:',...
        {' '},num2str(results{ii}.Cost_lambda),{' '},...
        '\cdot ||(vf,thk)||^2'));
    
    saveas(h,sprintf('FIG%d.png',ii));
    waitforbuttonpress;
end


