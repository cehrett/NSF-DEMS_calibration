function output = dual_calib_example_fn(...
    x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,ymean,ysd)
% This function is a univariate objective function with state-aware true
% optimum given by theta1 = 4/3 * (xx-1).

% Return inputs to original scale:
x = x * xrange + xmin    ; 
t1 = t1 * t1range + t1min ;
t2 = t2 * t2range + t2min ;

% Define function
f = @(x,t1,t2) x.*(t2.^(t1-1) .* exp(-0.75 * t2) + 1) .^ (-1) ;

% Get output on original scale
output_os = f(x(:),t1(:),t2(:)) ; 

% Standardize model outputs
output = (output_os - ymean) ./ ysd ;

end