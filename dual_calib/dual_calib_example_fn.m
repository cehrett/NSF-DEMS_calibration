function output = dual_calib_example_fn(...
    x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,ymean,ysd,discrep)
% This function is a univariate objective function with state-aware true
% optimum given by theta2 = 4/3 * (theta1-1).
% discrep is a nonneg integer identifying which discrepancy to use.

switch nargin
    case 3, rescale_inputs = false; discrep = 0;
        t1 = xmin; t2 = xrange; ymean = 0; ysd = 1;
    case 4, rescale_inputs = false; discrep = t1;
        t1 = xmin; t2 = xrange; ymean = 0; ysd = 1; 
        xmin = min(x); xrange=range(x);
    case 11, rescale_inputs = true; discrep = 0;
    case 12, rescale_inputs = true;
    otherwise, error('Incorrect number of inputs');
end

% Return inputs to original scale:
if rescale_inputs
    x = x * xrange + xmin    ; 
    t1 = t1 * t1range + t1min ;
    t2 = t2 * t2range + t2min ;
end


% Define objective function
obj_fn = @(xx,tt1,tt2) xx.*(tt2.^(tt1-1) .* exp(-0.75 * tt2) + 1).^(-1);

% Combine obj_fn with discrepancy
switch discrep % Select which discrepancy to use
    case 0 % No discrepancy
        f = obj_fn;
    case 1 % Multiplicative discrepancy, depends only on x
        c=1.5; %Raising/lowering c makes discrepancy more/less aggressive
        mult_disc_fn =@(xx) 1 - c*(xx-xmin).*(xx-(xmin+xrange))./xx;
        f = @(xx,tt1,tt2) ...
            mult_disc_fn(xx) .* obj_fn(xx,tt1,tt2);
    case 2 % Multiplicative discrepancy, depends only on x, same as 1 but c
        c=3.5; %Raising/lowering c makes discrepancy more/less aggressive
        mult_disc_fn =@(xx) 1 - c*(xx-xmin).*(xx-(xmin+xrange))./xx;
        f = @(xx,tt1,tt2) ...
            mult_disc_fn(xx) .* obj_fn(xx,tt1,tt2);
    case 3 % When theta1=2, this discrepancy does not change optimal theta2
        c = .15; %Raising/lowering c makes discrepancy more/less aggressive
        d=0.075; %Raising/lowering c makes discrepancy more/less aggressive
        add_disc_fn =@(xx,tt2) -c.*(xx-.5).*(xx-1).*(tt2 - 4/3).^2+d;
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 4 % Same as 3 except for c
        c = .65; %Raising/lowering c makes discrepancy more/less aggressive
        d=0.075; %Raising/lowering c makes discrepancy more/less aggressive
        add_disc_fn =@(xx,tt2) -c.*(xx-.5).*(xx-1).*(tt2 - 4/3).^2+d;
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 5 % This additive discrepancy changes optimal theta2
        % The new opt. theta2 is (approx) 1, when theta1=2 and c=0.055
        c = 0.055 ; % This controls the size of the additive discrepancy
        d = 0; %Raising/lowering c makes discrepancy more/less aggressive
        add_disc_fn = @(xx,tt2) c * xx .* tt2 + d; % Works better w/ d=0.1
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 6 % Same as 5 except for d
        % The new opt. theta2 is (approx) 1, when theta1=2 and c=0.055
        c = 0.055 ; % This controls the size of the additive discrepancy
        d = 0.1; %Raising/lowering d makes discrepancy more/less aggressive
        add_disc_fn = @(xx,tt2) c * xx .* tt2 + d; % Works better w/ d=0.1
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
end        

% Get output on original scale
output_os = f(x,t1,t2) ; 

% Standardize model outputs
output = (output_os - ymean) ./ ysd ;

end