% Get the covariance matrix for two arrays input and input_s

function R = gp_cov(rho,input,input_s,lambda,verbose)

% Initialize the covariance matrix
R = ones(size(input,1),size(input_s,1));

% If either input is empty, then output an empty cov matrix
if min(size(input))==0 || min(size(input_s))==0
    R=0*R;
    return;
end

% If lambda is Inf, then skip the pointless work and output 0 matrix
if lambda == Inf
    R = 0*R;
    return;
end

msg=0; % used for verbose output
% Now loop through all the elements of omega.
for ii = 1:numel(rho)
    inputmat = repmat(input(:,ii),1,size(input_s,1));
    inputmat_s = repmat(input_s(:,ii)',size(input,1),1);
    Exp_input = 4 * (inputmat-inputmat_s).^2;
    rhomat = rho(ii) * ones(size(R)); % Get a matrix of the current omega
    R = R .* rhomat.^Exp_input ; 
    if verbose
        percdone = 100*ii/numel(rho);
        fprintf(repmat('\b',1,msg));
        msg=fprintf('%3.2f%% done\n',percdone);
    end
end

% Use marginal precision:
R = R/lambda;

end
