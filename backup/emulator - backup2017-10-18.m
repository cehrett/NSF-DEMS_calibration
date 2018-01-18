function em = emulator(sim_dat,pred_pts,omega,rho,lambda,verbose)

dat=sim_dat;

% Rescale observered simulator inputs
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

% Rescale prediction inputs
temp_pred  = pred_pts(:,1);
VF_pred    = pred_pts(:,2);
thick_pred = pred_pts(:,3);
temp_pred_min   = min(temp_pred);
temp_pred_range = range(temp_pred);
VF_pred_min   = min(VF_pred);
VF_pred_range = range(VF_pred);
thick_pred_min   = min(thick_pred);
thick_pred_range = range(thick_pred);

temp_pred01 = (temp_pred-min(temp_pred))./range(temp_pred);
VF_pred01 = (VF_pred-min(VF_pred))./range(VF_pred);
thick_pred01 = (thick_pred-min(thick_pred))./range(thick_pred);


% Standardize observed simulator outputs
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

% Package things up for emulation
sim_des     = [temp_sim01 VF_sim01 thick_sim01];
sim_res     = [defl_sim_std rot_sim_std cost_sim_std];
pred_des    = [temp_pred01 VF_pred01 thick_pred01];

% Prepare prediction points and simulator design points
ones_vec = ones(size(pred_pts,1),1);
xp = [ 0*ones_vec pred_des ; .5 * ones_vec pred_des ; 
    ones_vec pred_des ];
ones_vec = ones(size(sim_des,1),1);
x  = [ 0*ones_vec sim_des; .5* ones_vec sim_des ; 
    ones_vec sim_des ];

% This will store the predicted response
ypred = zeros(size(pred_pts,1),1);

% Get Cov matrices
if verbose==true
    fprintf('Getting three covariance matrices... ')
end
Sigma_xx  = gp_cov(omega,x(:,1:2),x(:,1:2),...
    rho,x(:,3:end),x(:,3:end),lambda);
if verbose==true
    fprintf('1,... ')
end
Sigma_xpx  = gp_cov(omega,xp(:,1:2),x(:,1:2),...
    rho,xp(:,3:end),x(:,3:end),lambda);
if verbose==true
    fprintf('2,... ')
end
Sigma_xxp  = Sigma_xpx';
Sigma_xpxp  = gp_cov(omega,xp(:,1:2),xp(:,1:2),...
    rho,xp(:,3:end),xp(:,3:end),lambda);
if verbose==true 
    fprintf('3.\n\n')
end
WN = eye(size(Sigma_xx)) * 10^(-6);
% Now get the mean and covariance of the gp.
inv_Sig_xx = inv(Sigma_xx + WN);
rc=rcond(inv_Sig_xx);
%fprintf('rcond(inv(Sigma_xx+WN))=%.3g\n\n',rc);

% Get mean and covariance at prediction points
mu = Sigma_xpx * inv_Sig_xx * sim_res(:);
cov_gp = Sigma_xpxp - Sigma_xpx * inv_Sig_xx * Sigma_xxp;

% Pack up and leave
em.Sigma_xx = Sigma_xx;
em.mu       = mu;
em.cov_gp   = cov_gp;
em.rc       = rc;
em.xp       = xp;
em.sim_des  = sim_des;
em.sim_res  = sim_res;
em.orig_dat = dat;

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
