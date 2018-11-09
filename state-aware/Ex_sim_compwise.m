function output = Ex_sim_compwise(xx,theta1,theta2,...
    input_cntrl_min,input_cntrl_range,input_calib_min,input_calib_range,...
    output_mean,output_sd,which_sa)
% This function organizes the inputs xx,theta1,theta2 for evaluation of the
% function Ex_sim.m. This function also handles the rescaling of the model
% inputs, to return them to their original scale. The output is given as a
% vector, rather than as a three-column matrix (as is the direct output of
% Ex_sim.m)

% Get distinct control input values
c_vals = unique(xx(:,3),'stable');

switch mat2str(which_sa)
    case '[0 0]'
        t1 = repmat(theta2(1),length(c_vals),1) ;
        t2 = repmat(theta2(2),length(c_vals),1) ;
    case '[1 0]'
        t1 = theta1(:) ; 
        t2 = repmat(theta2,length(t1),1) ; 
    case '[0 1]' 
        t2 = theta1(:) ;
        t1 = repmat(theta2,length(t2),1) ; 
end

% Return inputs to original scale:
c_vals = c_vals * input_cntrl_range + input_cntrl_min    ; 
t1 = t1 * input_calib_range(1) + input_calib_min(1) ;
t2 = t2 * input_calib_range(2) + input_calib_min(2) ;

% Get model outputs (on original scale)
ct = [ c_vals t1 t2 ] ;
mat_out = Ex_sim(ct);

% Standardize model outputs
mat_out = (mat_out - output_mean) ./ output_sd ;

output = mat_out(:);

end