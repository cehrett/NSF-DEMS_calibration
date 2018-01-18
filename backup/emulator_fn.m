function em = emulator_fn(sim_dat,pred_pts,omega,rho)

dat=sim_dat;

% Rescale inputs
temp_sim  = dat(:,1);
VF_sim    = dat(:,2);
thick_sim = dat(:,3);
temp_sim_min   = min(temp_sim);
temp_sim_range = range(temp_sim);
VF_sim_min   = min(VF_sim);
VF_sim_range = range(VF_sim);
thick_sim_min   = min(thick_sim);
thick_sim_range = range(thick_sim);

temp_sim01 = (temp_sim-min(temp_sim))./range(temp_sim);
VF_sim01 = (VF_sim-min(VF_sim))./range(VF_sim);
thick_sim01 = (thick_sim-min(thick_sim))./range(thick_sim);

% Standardize outputs
defl_sim = dat(:,4);
rot_sim  = dat(:,5);
cost_sim = dat(:,6);
defl_sim_mean = mean(defl_sim);
defl_sim_sd   = std(defl_sim);
rot_sim_mean = mean(rot_sim);
rot_sim_sd   = std(rot_sim);
cost_sim_mean = mean(cost_sim);
cost_sim_sd   = std(cost_sim);

defl_sim_std = (defl_sim - defl_sim_mean)/defl_sim_sd;
rot_sim_std  = (rot_sim - rot_sim_mean)/rot_sim_sd;
cost_sim_std = (cost_sim - cost_sim_mean)/cost_sim_sd;

sim_des     = [temp_sim01 VF_sim01 thick_sim01];
sim_res     = [defl_sim_std rot_sim_std cost_sim_std];


% Prepare prediction points
ones_vec = ones(size(pred_pts,1),1);
pred_pts = [ 0*ones_vec pred_pts ; .5 * ones_vec pred_pts ; 
    ones_vec pred_pts ];

% This will store the predicted response
ypred = zeros(size(pred_pts,1),1);

xp = pred_pts;
fprintf('Getting four covariance matrices... ')
Sigma_xx  = gp_cov(omega,x(:,1:2),x(:,1:2),...
    rho,x(:,3:end),x(:,3:end));
fprintf('1,... ')
Sigma_xpx  = gp_cov(omega,xp(:,1:2),x(:,1:2),...
    rho,xp(:,3:end),x(:,3:end));
fprintf('2,... ')
Sigma_xxp  = gp_cov(omega,x(:,1:2),xp(:,1:2),...
    rho,x(:,3:end),xp(:,3:end));
fprintf('3.\n\n')
Sigma_xpxp  = gp_cov(omega,pred_pts(:,1:2),pred_pts(:,1:2),...
    rho,pred_pts(:,3:end),pred_pts(:,3:end));
WN = eye(size(Sigma_xx)) * 10^(-6);
% Now get the mean and covariance of the gp.
inv_Sig_xx = inv(Sigma_xx + WN);
rc=rcond(inv_Sig_xx);
fprintf('rcond(inv(Sigma_xx+WN))=%.3g',rc);
mu = Sigma_xpx * inv_Sig_xx * sim_res(:);
cov_gp = Sigma_xpxp - Sigma_xpx * inv_Sig_xx * Sigma_xxp;

% Pack up and leave
em.Sigma_xx = Sigma_xx;
em.mu       = mu;
em.cov_gp   = cov_gp;
em.rc       = rc;

stds.temp_sim_min    = temp_sim_min;
stds.temp_sim_range  = temp_sim_range;
stds.VF_sim_min      = VF_sim_min;
stds.VF_sim_range    = VF_sim_range;
stds.thick_sim_min   = thick_sim_min;
stds.thick_sim_range = thick_sim_range;

stds.defl_sim_mean = defl_sim_mean;
stds.defl_sim_sd   = defl_sim_sd;
stds.rot_sim_mean  = rot_sim_mean;
stds.rot_sim_sd    = rot_sim_sd;
stds.cost_sim_mean = cost_sim_mean;
stds.cost_sim_sd   = cost_sim_sd;

em.stds = stds;

end
