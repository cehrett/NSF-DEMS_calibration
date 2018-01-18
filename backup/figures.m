% Figures from emulator; this code can be run only after emulator.m

% Unpack em struct
mu           = em.mu;
cov_gp       = em.cov_gp;
xp           = em.xp;
% sim_des      = em.sim_des;
% sim_res      = em.sim_res;
% sim_dat_orig = em.orig_dat;

% Return inputs to original scale
temp_sim_min    = tdat.input_mins(2);
temp_sim_range  = tdat.input_ranges(2);
VF_sim_min      = tdat.input_mins(3);
VF_sim_range    = tdat.input_ranges(3);
thick_sim_min   = tdat.input_mins(4);
thick_sim_range = tdat.input_ranges(4);

xp_orig = xp; % This will become original scale
xp_orig(:,2) = xp_orig(:,2)*temp_sim_range + temp_sim_min;
xp_orig(:,3) = xp_orig(:,3)*VF_sim_range + VF_sim_min;
xp_orig(:,4) = xp_orig(:,4)*thick_sim_range + thick_sim_min;
% ones_vec = ones(size(xp,1),1);
% xp_orig = [ 0 * ones_vec xp_orig ; .5 * ones_vec xp_orig ; ...
%     1 * ones_vec xp_orig];

temp_pred  = unique(xp_orig(:,2));
VF_pred    = unique(xp_orig(:,3));
thick_pred = unique(xp_orig(:,4));

% Get values to put emulator output back on original scale
defl_mean = tdat.output_means(1,:)';
defl_sd   = tdat.output_sds(1,:)';
rot_mean  = tdat.output_means(2,:)';
rot_sd    = tdat.output_sds(2,:)';
cost_mean = tdat.output_means(3,:)';
cost_sd   = tdat.output_sds(3,:)';

%% Bring mu back to original scale
% defl_mean_s will store a vector the same length as mu in which each
% element is the mean at the temperature of the corresponding point in mu.
% defl_sd_s is similar. We repeat for rot and cost.
defl_mean_s = repelem(...
    defl_mean,length(mu)/length(temp_pred)/3);
defl_sd_s = repelem(...
    defl_sd,length(mu)/length(temp_pred)/3);
rot_mean_s = repelem(...
    rot_mean,length(mu)/length(temp_pred)/3);
rot_sd_s = repelem(...
    rot_sd,length(mu)/length(temp_pred)/3);
cost_mean_s = repelem(...
    cost_mean,length(mu)/length(temp_pred)/3);
cost_sd_s = repelem(...
    cost_sd,length(mu)/length(temp_pred)/3);
means = [ defl_mean_s ; rot_mean_s ; cost_mean_s ] ;
sds   = [ defl_sd_s ; rot_sd_s ; cost_sd_s ] ;
mu_orig = mu .* sds + means;

% Let's look at the deflection surface for just temp and VF, 
% at one thickness (10).
figure();
desired_thickness = 10;
idx1 = xp_orig(:,1)==0 & xp_orig(:,4)==desired_thickness;
[EE,GG] = meshgrid(temp_pred,VF_pred);
mu1 = reshape(mu_orig(idx1),ms,length(temp_pred));
surf(EE,GG,mu1);
xlabel('Temperature');
ylabel('VF');
zlabel('Tip defl.');
title(['Predictions and FE model observations at thickness = ' ...
    num2str(desired_thickness) ]);
%Let's compare the surface to observations at this thickness.
idxc1 = raw_dat(:,3)==desired_thickness;
hold on;
plot3(raw_dat(idxc1,1),raw_dat(idxc1,2),...
    raw_dat(idxc1,4),'r.','MarkerSize', 20);
hold off;

% Now again for twist angle.
figure();
idx2 = xp(:,1)==.5 & xp(:,4)==xp(2,4);
mu2 = reshape(mu(idx2),ms,ms);
surf(EE,GG,mu2);
xlabel('Temperature');
ylabel('VF');
zlabel('Twist angle');
%Let's compare the surface to observations at nearby thicknesses.
idxc2 = dat(:,3)==15;
hold on;
plot3(sim_des(idxc2,1),sim_des(idxc2,2),sim_res(idxc2,2),'ko');
hold off;