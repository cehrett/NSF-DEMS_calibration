function output = dual_calib_example_fn(...
    x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,ymean,ysd,discrep,...
    rescale_inputs,c,d)
% This function is a univariate objective function with state-aware true
% optimum given by theta2 = 4/3 * (theta1-1).
% discrep is a nonneg integer identifying which discrepancy to use.

if xmin == 0, error('Error: Need minima and ranges for all inputs'); end

switch nargin
    case 3, error('Error: Need minima and ranges for all inputs');
%         rescale_inputs = false; discrep = 0;
%         t1 = xmin; t2 = xrange; ymean = 0; ysd = 1; setcd = false;
    case 4, error('Error: Need minima and ranges for all inputs');
%         rescale_inputs = false; discrep = t1;
%         t1 = xmin; t2 = xrange; ymean = 0; ysd = 1; 
%         xmin = min(x); xrange=range(x); setcd = false;
    case 12, discrep = 0; setcd = false;
    case 13, setcd = false;
    case 15, setcd = true; disp('ye');
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
        if setcd == false % If c,d were not specified by the user:
            c=1.5; %Raise/lower c makes discrepancy more/less aggressive
        end
        mult_disc_fn =@(xx) 1 - c*(xx-xmin).*(xx-(xmin+xrange))./xx;
        f = @(xx,tt1,tt2) ...
            mult_disc_fn(xx) .* obj_fn(xx,tt1,tt2);
    case 2 % Multiplicative discrepancy, depends only on x, same as 1 but c
        if setcd == false % If c,d were not specified by the user:
            c=3.5; %Raising/lower c makes discrepancy more/less aggressive
        end
        mult_disc_fn =@(xx) 1 - c*(xx-xmin).*(xx-(xmin+xrange))./xx;
        f = @(xx,tt1,tt2) ...
            mult_disc_fn(xx) .* obj_fn(xx,tt1,tt2);
    case 3 % When theta1=2, this discrepancy does not change optimal theta2
        if setcd == false % If c,d were not specified by the user:
            c = .15; %Raise/lower c makes discrepancy more/less aggressive
            d=0.075; %Raise/lower d makes discrepancy more/less aggressive
        end
        add_disc_fn =@(xx,tt2) -c.*(xx-.5).*(xx-1).*(tt2 - 4/3).^2+d;
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 4 % Same as 3 except for c
        if setcd == false % If c,d were not specified by the user:
            c = .65; %Raise/lower c makes discrepancy more/less aggressive
            d=0.075; %Raise/lower d makes discrepancy more/less aggressive
        end
        add_disc_fn =@(xx,tt2) -c.*(xx-.5).*(xx-1).*(tt2 - 4/3).^2+d;
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 5 % This additive discrepancy changes optimal theta2
        % The new opt. theta2 is (approx) 1, when theta1=2 and c=0.055
        if setcd == false % If c,d were not specified by the user:
            c=0.055; %Raise/lower c makes discrepancy more/less aggressive
            d = 0;   %Raise/lower d makes discrepancy more/less aggressive
        end
        add_disc_fn = @(xx,tt2) c * xx .* tt2 + d; % Works better w/ d=0.1
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
    case 6 % Same as 5 except for d
        % The new opt. theta2 is (approx) 1, when theta1=2 and c=0.055
        if setcd == false % If c,d were not specified by the user:
            c=0.055; %Raise/lower c makes discrepancy more/less aggressive
            d = 0.1; %Raise/lower d makes discrepancy more/less aggressive
        end
        add_disc_fn = @(xx,tt2) c * xx .* tt2 + d; % Works better w/ d=0.1
        f = @(xx,tt1,tt2) obj_fn(xx,tt1,tt2) + add_disc_fn(xx,tt2);
end        

% Get output on original scale
output_os = f(x,t1,t2) ; 

% Standardize model outputs
output = (output_os - ymean) ./ ysd ;

end