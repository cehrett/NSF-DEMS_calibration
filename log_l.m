% Log likelihood function for data

function L = log_l (Eta,X,Theta,omega,rho)

R = gp_cov(omega,X,X,rho,Theta,Theta);

L = -1/2 * log(det(R)) - lambda/2 * (Eta - mu)' * inv(R) * (Eta-mu)

end