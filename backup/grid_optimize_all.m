%Grid_optimize

function res = grid_optimize_all(dat, omega1, omega2, rho1, rho2, lambda)
% Data, length of grid, parameters.

% 
count = 0 ; % To keep track of progress
msg=0; % Will be used for progress indicator

% Separate data into training and validation sets.
idx = zeros(size(dat,1),1);
for ii = 1:length(idx)
    idx(ii) = mod(ii,10)==0;
end
idx_tr = idx==0;
idx    = idx==1; %converts to logical
dat_tr  = dat(idx_tr,:);
dat_val = dat(idx,:);

best_MSPE=10^10; % High value, replace with lowest at each step

for ii=1:length(omega1)
    for jj=1:length(omega2)
        for kk=1:length(rho1)
            for ll=1:length(rho2)
                for mm=1:length(lambda)
                    omega  = [omega1(ii) omega2(jj)];
                    rho    = [rho1(kk) rho2(ll)];
                    lam = lambda(mm);
                    em = emulator(dat_tr,dat_val,omega,rho,lam,false);

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
                    
                    % Print update and, if applicable, record new bests
                    count = count + 1; pp = 100*count / (length(omega1)*length(omega2)*length(rho1)*length(rho2)*length(lambda));
                    fprintf(repmat('\b',1,msg));
                    if MSPE < best_MSPE
                        fprintf(repmat('\b',1,msg));
                        fprintf('Best MSPE: %8.4g\n',MSPE);
                        fprintf('Best rho: %6.3f,%6.3f\n',rho(1),rho(2));
                        fprintf('Best omega: %6.3f,%6.3f\n',omega(1),omega(2));
                        fprintf('Best lambda: %6.4f\n',lam);
                        fprintf('Recorded at %3.1f%% completion\n\n\n',pp);
                        msg=0;
                        res.omega  = omega;
                        res.rho    = rho;
                        res.lambda = lam;
                        res.mspe   = MSPE;
                        best_MSPE  = MSPE;
                    end
                                        
                    
                    msg = fprintf('Last MSPE: %8.4g\n',MSPE);
                    msg = msg + fprintf('Last rho: %6.3f,%6.3f\n',rho(1),rho(2));
                    msg = msg + fprintf('Last omega: %6.3f,%6.3f\n',omega(1),omega(2));
                    msg = msg + fprintf('Last lambda: %6.4f\n\n\n',lam);
                    msg = msg + fprintf('Grid optimization routine %3.1f%% complete\n\n',pp);
                end
            end
        end
        
    end
end

end