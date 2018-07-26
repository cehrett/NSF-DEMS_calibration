function h = calib_heatmap (des_obs_os,post_theta,mfa,mea,ms)
% Generate a heatmap of proximity to desired observation using direct data

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 
if ~exist('post_theta','var')
    post_theta=0
end
if ~exist('mfa','var')
    mfa=0.5;
end
if ~exist('mea','var')
    mea=0.5;
end
if ~exist('ms','var')
    ms=5;
end

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);

% Load MCMC samples
load([dpath,'Example\Ex_results\'...
    '2018-07-11_discrepancy_true_fn_set_lambda_delta_1'],...
    'results');
% load([dpath,'Example\Ex_results\'...
%     '2018-07-12_discrepancy_true_fn_set_lambda_delta_1-256'],...
%     'results');

% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get true samples with output closest to 0 (Euclidean distance on
% standardized scale
% First put data on standardized scale
meanout = mean(results.settings.output_means');
sdout   = mean(results.settings.output_sds'  );
cost_std = (ctheta_output(:,6) - meanout(3))/...
    sdout(3);
defl_std = (ctheta_output(:,4) - meanout(1))/...
    sdout(1);
rotn_std = (ctheta_output(:,5) - meanout(2))/...
    sdout(2);

dd_outputs_std = [defl_std rotn_std cost_std];

% Get desired_obs on standardized scale
des_obs = (des_obs_os-meanout)./...
    sdout;

% Now get Euclidean norms of each standardized output
dd_dists = sqrt ( sum ( (dd_outputs_std-des_obs).^2 , 2 ) ) ;

% Now get MCMC sample output and put on standardized scale
outs = results.model_output.by_sample_true(results.settings.burn_in:end,:);
% Put on standardized scale:
cost_std = (outs(:,3) - meanout(3))/...
    sdout(3);
defl_std = (outs(:,1) - meanout(1))/...
    sdout(1);
rotn_std = (outs(:,2) - meanout(2))/...
    sdout(2);

mcmc_outputs_std = [ defl_std rotn_std cost_std ] ;

% Now get Euclidean norms of each standardized output
mcmc_dists = sqrt ( sum ( (mcmc_outputs_std-des_obs).^2 , 2 ) ) ;
 
% Take a look
% figure();
% scatter(linspace(1,length(mcmc_dists),length(dd_dists)),dd_dists);
% hold on;
% scatter(1:length(mcmc_dists),mcmc_dists);
% 
% % Now take a 3d look at all outputs versus the close direct data outputs
cutoff = quantile(mcmc_dists,.95); % cutoff for close dd output
close_dd_idx = dd_dists <= cutoff; % index of close dd outputs
close_dd_outputs = ctheta_output(close_dd_idx,4:6) ; % close dd outputs
% figure();
% scatter3(outs(:,1),outs(:,2),outs(:,3),20); axis vis3d; hold on;
% scatter3(...
%     close_dd_outputs(:,1),close_dd_outputs(:,2),close_dd_outputs(:,3));

% Now take a look at all calib settings at mcmc outputs vs close dd outputs

close_dd_theta = ctheta_output(close_dd_idx,2:3);
% figure(); scatterhist(samps(:,1),samps(:,2));
% figure(); scatterhist(close_dd_theta(:,1),close_dd_theta(:,2));

% Now get a scatterhist of mcmc theta draws with, behind it, all direct
% data theta values colored by Euclidean distance of the standardized
% output to the zero point.
h=figure(); colormap(flipud(jet));
if ~isequal(post_theta,0) % ie if posterior draws were supplied
    samps = post_theta;
    sh=scatterhist(samps(:,1),samps(:,2),'Marker','.'); 
    hold on; 
end
xlim([0 3]); ylim([0 6]);
scatter(ctheta_output(:,2),ctheta_output(:,3),2,dd_dists); hold on;
colorbar('East');
if ~isequal(post_theta,0) % ie if posterior draws were supplied
    scatter(samps(:,1),samps(:,2),ms,'.g','MarkerFaceAlpha',mfa,...
    'MarkerEdgeAlpha',mea);
end
end