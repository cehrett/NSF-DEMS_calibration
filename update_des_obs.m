function [new_des_obs,new_lambda_delta] = ...
    update_des_obs (results,desired_obs)
% Given results of calibration, this function provides settings for a new
% round of calibration with an updated desired observation, and
% corresponding updated estimate of lambda_delta, the marginal precision of
% the discrepancy function.
% 
% The function works by locating the closest 100 samples from the results,
% performing regression to get a hyperplane, then updating the desired
% observation to lie in the line segment connecting the original desired
% observation and this hyperplane.
%
% The distance of the new desired observation is half that of the old one.

des_obs_os = desired_obs ; %results.settings.desired_obs_orig ; 
des_obs = (des_obs_os - mean(results.settings.output_means,2)')./...
    mean(results.settings.output_sds,2)';
n=size(results.samples,1);
if isfield(results.model_output,'by_sample_true')
    outs_os = results.model_output.by_sample_true;
else
    outs_os = results.model_output.by_sample_est ; 
end
outs = (outs_os - mean(results.settings.output_means,2)')./...
    mean(results.settings.output_sds,2)';

%% Get distances of model outputs from desired observation, find close ones
dists = sum( (outs - des_obs).^2,2 );
close_idx = dists < quantile(dists,100/n) ;
close_outs = outs(close_idx,:);
mean_close_out = mean(close_outs);

%% Get new desired observation direction and distance
dirvec_unnormed = mean_close_out - des_obs;
dist_to_des_obs = norm(dirvec_unnormed);
dirvec = dirvec_unnormed / norm(dirvec_unnormed);
dists_close_outs = sqrt(sum( (mean_close_out - close_outs).^2, 2));
new_dist = 3 * std(dists_close_outs);

%% Get new desired observation and distance
new_des_obs_ss = mean_close_out - new_dist * dirvec ;
% Convert back to original scale
new_des_obs = new_des_obs_ss .* mean(results.settings.output_sds,2)' + ...
    mean(results.settings.output_means,2)' ;
new_lambda_delta = 1/new_dist^2;

% %scratch
% close_outs_dists = sqrt(sum( (close_outs - des_obs).^2,2));
% scatter3(close_outs(:,1),close_outs(:,2),close_outs(:,3));
% hold on;
% scatter3(des_obs(:,1),des_obs(:,2),des_obs(:,3));
% line([des_obs(1) mean_close_out(1)],[des_obs(2) mean_close_out(2)],...
%     [des_obs(3) mean_close_out(3)]);
% scatter3(new_des_obs(:,1),new_des_obs(:,2),new_des_obs(:,3));
% scatter3(do(:,1),do(:,2),do(:,3));


end