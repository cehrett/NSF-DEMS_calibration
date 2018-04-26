function em = emulator(sim_dat_input,sim_dat_output,...
    pred_pts,omega,rho,lambda,Sigma_xx,inv_Sig_xx,verbose)
% This function outputs a covariance matrix for the input control and
% calibration vectors.

% sim_dat_input should be rescaled to [0,1] and sim_dat_output should be
% standardized. pred_pts on the same scales. Both should already
% have a dummy variable corresponding to output index as their first
% column.
% The covariance matrix of the training points may optionally be supplied.
% if not available, set argument Sigma_xx to be 0.

% Package things up for emulation
x           = sim_dat_input;
sim_res     = sim_dat_output;
xp          = pred_pts;

% This will store the predicted response
ypred = zeros(size(pred_pts,1),1);

% Get Cov matrices
WN_xx = eye(size(x,1)) * 10^(-2); % For numerical stabilization
WN_xpxp = eye(size(pred_pts,1)) * 10^(-2); % For numerical stabilization
if verbose==true
    fprintf('Getting three covariance matrices... ')
end

% Get Cov matrix of training points only if it is not supplied
if Sigma_xx == 0
    Sigma_xx  = gp_cov(omega,x(:,1:size(x,2)-2),x(:,1:size(x,2)-2),...
        rho,x(:,size(x,2)-1:end),x(:,size(x,2)-1:end),lambda,verbose) +...
        WN_xx;
    inv_Sig_xx = inv(Sigma_xx);
    rc=rcond(inv_Sig_xx);
end

if verbose==true
    fprintf('1,... ')
end
Sigma_xpx  = gp_cov(omega,xp(:,1:size(x,2)-2),x(:,1:size(x,2)-2),...
    rho,xp(:,size(x,2)-1:end),x(:,size(x,2)-1:end),lambda,verbose);
if verbose==true
    fprintf('2,... ')
end
Sigma_xxp  = Sigma_xpx';
Sigma_xpxp  = gp_cov(omega,xp(:,1:size(x,2)-2),xp(:,1:size(x,2)-2),...
    rho,xp(:,size(x,2)-1:end),xp(:,size(x,2)-1:end),lambda,verbose) +...
    WN_xpxp;
if verbose==true 
    fprintf('3.\n\n')
end

%fprintf('rcond(inv(Sigma_xx+WN))=%.3g\n\n',rc);

% Get mean and covariance at prediction points
mu = Sigma_xpx * inv_Sig_xx * sim_res(:);
cov_gp = Sigma_xpxp - Sigma_xpx * inv_Sig_xx * Sigma_xxp;

% Pack up and leave
em.Sigma_xx = Sigma_xx;
em.mu       = mu;
em.cov_gp   = cov_gp;
% em.rc       = rc;
em.xp       = xp;
% em.sim_des  = sim_des;
% em.sim_res  = sim_res;
% em.orig_dat = dat;

%%% Commenting the below out because it is now done elsewhere
% stds.temp_sim_min    = temp_sim_min;
% stds.temp_sim_range  = temp_sim_range;
% stds.VF_sim_min      = VF_sim_min;
% stds.VF_sim_range    = VF_sim_range;
% stds.thick_sim_min   = thick_sim_min;
% stds.thick_sim_range = thick_sim_range;
% 
% stds.defl_sim_mean = defl_sim_mean;
% stds.defl_sim_sd   = defl_sim_sd;
% stds.rot_sim_mean  = rot_sim_mean;
% stds.rot_sim_sd    = rot_sim_sd;
% stds.cost_sim_mean = cost_sim_mean;
% stds.cost_sim_sd   = cost_sim_sd;
% 
% em.stds = stds;

end
