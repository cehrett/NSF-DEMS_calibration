function qvals = nearquant(q,x,obs,locs,k)
% This function calculates quantile q at each point of x using k nearest
% neighbors observations (obs) at each observed location (locs).

% Use k-nearest neighbors
whichlocs = knnsearch(locs,x','K',k);

% Use all within distance d
%whichlocs = reshape(all([

nearobs   = obs(whichlocs);
qvals     = quantile(nearobs',q);

end