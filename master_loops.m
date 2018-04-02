% Master file for looped runs

clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');


%% User defined values
M = 1e3;
desired_obs = [.65 0.077 96];
desired_obs = [.65 96];
desired_obs = [.65 0.77];
which_outputs = [ 1 1 0 ] ; %Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Joint prop for theta, joint prop for obs var, prior on obs var
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);


post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
    settings.sim_xt,settings.eta,settings.output_sds,...
    settings.output_means,settings.omega,settings.rho,...
    settings.lambda,which_outputs);
results = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
    'init',samples(1,:),'desired_data',desired_obs,...
    'sigma2_prior',settings.log_sigma2_prior,...
    'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
    'proposal',settings.proposal,'nugsize',settings.nugsize,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
    'post_mean_out',post_mean_out);

%save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
%'results_1d_homosked2sdposdObsPrior'],...
%'results');

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
'results_0d_heteroskedPriorObs_hasBndCnds'],...
'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parallel processing loop
n=10;
design = lhsdesign(n,2) .* [ .74 224 ] ;
p = 3 ; % Number of chains at each desired obs. to run in parallel
% Cell to store data
results_par = cell(p,1);
parpool(p);
M=1e4;
all_results2 = cell(n,1);
% Begin from here if interrupted
load('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\30_MCMCs');

% Note: the range of the design was chosen by selecting 0 as the minimum
% values and for the maximum values taking the midpoint of the plausible
% ranges supplied by Evan for deflection and cost.

for jj=6:10

    parfor ii=1:p
        desired_obs = design(jj,:);
        settings = MCMC_settings(M,desired_obs);
        settings.doplot = false;
        [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
        post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
            settings.sim_xt,settings.eta,settings.output_sds,...
            settings.output_means,settings.omega,settings.rho,...
            settings.lambda);
        results = struct('samples',samples,'sigma2',sigma2_rec,...
            'Sigma',Sigma,'init',settings.init_theta,...
            'desired_data',desired_obs,...
            'sigma2_prior',settings.log_sigma2_prior,...
            'omega_rho_lambda',...
            [settings.omega settings.rho settings.lambda],...
            'proposal',settings.proposal,'nugsize',settings.nugsize,...
            'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
            'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
            'post_mean_out',post_mean_out);
        results_par{ii} = results;
    end
    
    all_results2{jj} = results_par;
    
    fprintf('COMPLETED LOOP %g/%g',jj,n)
    
    save('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\30_MCMCs');
    
end

save('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\30_MCMCs');

% Take a look at results
logit = @(x) log(x./(1-x));
h=figure('rend','painters','pos',[10 10 1200 800]);
for jj=1:n
    for ii=1:p
        subplot(2,3,1);
        plot(all_results{jj}{ii}.samples(burn_in:end,1),'ko');
        title('Volume fraction');
        subplot(2,3,2);
        plot(all_results{jj}{ii}.samples(burn_in:end,2),'ko');
        title('Thickness');
        subplot(2,3,4);
        plot(all_results{jj}{ii}.sigma2(burn_in:end,1),'ko');
        title('Deflection \sigma^2');
        subplot(2,3,5);
        plot(all_results{jj}{ii}.sigma2(burn_in:end,2),'ko');
        title('Cost \sigma^2');
        
        subplot(2,3,6);
        lbsamps = logit(all_results{jj}{ii}.samples(burn_in:end,:));
        plot(lbsamps(:,1),lbsamps(:,2),'ko');
        hold on;
        rr = mvnrnd(mean(lbsamps),all_results{jj}{ii}.Sigma,150);
        plot(rr(:,1),rr(:,2),'r.');
        hold off;
        title('(Logit) samples with proposal cov.');
        
        suptitle(strcat('Design point',{' '},num2str(jj),', Chain',...
            {' '},num2str(ii)));
        
        hh = subplot(2,3,3);
        vfacf = acf(all_results{jj}{ii}.samples(burn_in:end,1),50);
        thkacf = acf(all_results{jj}{ii}.samples(burn_in:end,2),50);
        vfacf = vfacf(50); thkacf=thkacf(50);
        cla;
        title(' ')
        line1Str = strcat(['VF posterior mean:' ' ' num2str(...
            mean(all_results{jj}{ii}.samples(burn_in:end,1)))]);
        line2Str = strcat(['Thk posterior mean:' ' ' num2str(...
            mean(all_results{jj}{ii}.samples(burn_in:end,2)))]);
        line3Str = strcat(['50 lag ACF for VF: ' ' ' num2str(vfacf)]);
        line4Str = strcat(['50 lag ACF for Thk:' ' ' num2str(thkacf)]);
        line5Str = strcat(['Desired deflection):' ' ' ...
            num2str(all_results{jj}{ii}.desired_data(1))]);
        line6Str = strcat(['Desired cost):' ' ' ...
            num2str(all_results{jj}{ii}.desired_data(2))]);
        xl = xlim(hh); 
        xPos = 0; 
        yl = ylim(hh); 
        yPos = yl(1) + diff(yl) / 2; 
        t = text(xPos, yPos, sprintf('%s\n%s\n\n%s\n%s\n\n%s\n%s', ...
            line1Str, line2Str, line3Str, line4Str,line5Str,line6Str), 'Parent',h);
        %set(t, 'HorizontalAlignment', 'center');
        %title('info')
        axis off
        waitforbuttonpress;
        %saveas(h,sprintf('FIG%d-%d.png',jj,ii));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Joint proposal for theta, prior put on obs variance
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior(M,burn_in,sim_xt,eta,...
    obs_x,y,sigma2,log_sigma2_prior,out_of_range,init_theta,omega,...
    rho,lambda,proposal,nugsize);

results = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
'init',init_theta,'desired_data',desired_obs,...
'sigma2_prior',log_sigma2_prior,...
'omega_rho_lambda',[omega rho lambda],'proposal',proposal,...
'nugsize',nugsize);

save(['.\NSF DEMS\Phase 1\'...
'results_z0_univObsPrior1'],...
'results');
