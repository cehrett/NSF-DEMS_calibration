function output = Univar_objective(xx,theta1,~,...
    input_cntrl_min,input_cntrl_range,input_calib_min,input_calib_range,...
    output_mean,output_sd,~)
% This function is a univariate objective function with state-aware true
% optimum given by theta1 = 4/3 * (xx-1).

% Return inputs to original scale:
xx = xx * input_cntrl_range + input_cntrl_min    ; 
t1 = theta1 * input_calib_range + input_calib_min ;

% Define function
f = @(x,t) (t.^(x-1) .* exp(-0.75 * t) + 1) .^ (-1) ;

% Get output on original scale
output_os = f(xx(:),t1(:)) ; 

% Standardize model outputs
output = (output_os - output_mean) ./ output_sd ;

end