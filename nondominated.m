function [nondoms,ndidx] = nondominated(results,tol)
% This function takes as input a set of outcomes and returns those that are
% nondominated, ie Pareto efficient. Minimization is assumed to be optimal.
% tol is an optional variable which loosens the idea of nondominance. It
% should be a vector of length equal to the number of columns in results. 

% tol is optional so set it to 0 if it isn't there
if ~exist('tol','var')
    tol = 0;
end

% This will hold the nondominated outcomes:
nondoms = [];

% This will hold the indices of those outcomes:
ndidx = [];

% Supplement results with indices:
results = [results [1:size(results,1)]'];

% The loop whittles down the results, keeping nondominated outcomes and
% throwing away outcomes that they dominate, until there are no results
% left.
while size(results,1) ~= 0 %

    % Get minimum on first objective
    [~,new_idx] = min(results(:,1));
    or_idx = results(new_idx,end); % Get idx of new min in original results

    nondoms = [nondoms ; results(new_idx,1:end-1) ] ; % Get new member 
                                                      % of nondoms
    ndidx   = [ndidx   ; or_idx                   ] ;

    % This tells us which outcomes are dominated by the lastest nondom elt.
    dominated_idx = ...
        all(results(:,1:end-1) - tol >= results(new_idx,1:end-1),2) ;
    % Add the index of the new nondom
    dominated_idx(new_idx) = true;

    % Throw away the dominated outcomes
    results = results(~dominated_idx,:);
    
    if rand < 0.0025
        fprintf('Still left: %d\n',size(results,1));
    end

end

end