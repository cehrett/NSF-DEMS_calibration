function output = dual_calib_example_fn(...
    x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,ymean,ysd,discrep)
% This function is a univariate objective function with state-aware true
% optimum given by theta1 = 4/3 * (xx-1).
% discrep is a boolean which, when true, adds a discrepancy to the output.

switch nargin
    case 11, discrep = false;
    case 12, discrep = boolean(discrep);
    otherwise, error('Incorrect number of inputs');
end

% Return inputs to original scale:
x = x * xrange + xmin    ; 
t1 = t1 * t1range + t1min ;
t2 = t2 * t2range + t2min ;

if discrep % Here's where we define the discrepancy
    disc_fn = @(xx,tt) ...
        xx + (xx - (xmin+xrange/2)).^2 .* tt/2 - tt/2 *xrange^2/4;
    mult = disc_fn(x,t2) ; 
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