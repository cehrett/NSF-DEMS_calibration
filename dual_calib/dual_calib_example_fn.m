function output = dual_calib_example_fn(...
    x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,ymean,ysd,discrep)
% This function is a univariate objective function with state-aware true
% optimum given by theta2 = 4/3 * (theta1-1).
% discrep is a boolean which, when true, adds a discrepancy to the output.

switch nargin
    case 3, rescale_inputs = false; discrep = false;
        t1 = xmin; t2 = xrange; ymean = 0; ysd = 1;
    case 4, rescale_inputs = false; discrep = boolean(t1);
        t1 = xmin; t2 = xrange; ymean = 0; ysd = 1; 
        xmin = min(x); xrange=range(x);
    case 11, rescale_inputs = true; discrep = false;
    case 12, rescale_inputs = true; discrep = boolean(discrep);
    otherwise, error('Incorrect number of inputs');
end

% Return inputs to original scale:
if rescale_inputs
    x = x * xrange + xmin    ; 
    t1 = t1 * t1range + t1min ;
    t2 = t2 * t2range + t2min ;
end

if discrep % Here's where we define the discrepancy
    disc_fn =@(xx) xx - 1.5*(xx-xmin).*(xx-(xmin+xrange));
    mult = disc_fn(x) ; 
else
    mult = x;
end

% Define objective function
f = @(x,t1,t2) mult.*(t2.^(t1-1) .* exp(-0.75 * t2) + 1) .^ (-1) ;

% Get output on original scale
output_os = f(x(:),t1(:),t2(:)) ; 

% Standardize model outputs
output = (output_os - ymean) ./ ysd ;

end