%Grid_optimize

function res = grid_optimize_omega(dat, omega1, omega2, rho1, rho2)
% Data, length of grid, parameters.

% 
count = 0 ; % To keep track of progress
rho = [rho1 rho2];

% Separate data into training and validation sets.
idx = zeros(size(dat,1),1);
for ii = 1:length(idx)
    idx(ii) = mod(ii,10)==0;
end
idx_tr = idx==0;
idx    = idx==1; %converts to logical
dat_tr  = dat(idx_tr,:);
dat_val = dat(idx,:);

for ii=1:length(omega1)
    for jj=1:length(omega2)
        omega = [omega1(ii) omega2(jj)];
        em = emulator(dat_tr,dat_val,omega,rho);
        
        % Rescale val point inputs
        temp_min    = em.stds.temp_sim_min;
        temp_range  = em.stds.temp_sim_range;
        VF_min      = em.stds.VF_sim_min;
        VF_range    = em.stds.VF_sim_range;
        thick_min   = em.stds.thick_sim_min;
        thick_range = em.stds.thick_sim_range;
        val_temp01  = (dat_val(:,1)-temp_min)/temp_range;
        val_VF01    = (dat_val(:,2)-VF_min)/VF_range;
        val_thick01 = (dat_val(:,3)-thick_min)/thick_range;
        val_inputs  = [ val_temp01 val_VF01 val_thick01 ];
        
        % Standardize val point outputs
        defl_mean   = em.stds.defl_sim_mean;
        defl_sd     = em.stds.defl_sim_sd;
        rot_mean    = em.stds.rot_sim_mean;
        rot_sd      = em.stds.rot_sim_sd;
        cost_mean   = em.stds.cost_sim_mean;
        cost_sd     = em.stds.cost_sim_sd;
        val_defl_std = (dat_val(:,4) - defl_mean)/defl_sd;
        val_rot_std  = (dat_val(:,5) - rot_mean)/rot_sd;
        val_cost_std = (dat_val(:,6) - cost_mean)/cost_sd;      
        val_outputs  = [val_defl_std val_rot_std val_cost_std];
        
        mu=em.mu;
        MSPE = sum((mu-val_outputs(:)).^2)/length(mu);
        
        res.omega1(ii,jj)        = omega(1);
        res.omega2(ii,jj)        = omega(2);
        res.omega_mspe(ii,jj)   = MSPE;
        count = count + 1; pp = 100*count / (length(omega1)*length(omega2));
        fprintf('Grid optimization routine %3.1f%% complete\n\n',pp);
        
    end
end

end