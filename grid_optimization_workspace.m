%Grid optimizer workspace

% 
clc; clear all; close all;

addpath('NSF DEMS\Phase 1\');

fprintf('Reading data from .xlsx...\n')
dat = xlsread("fe_results.xlsx");

%% Optimize rho
% Determine grid on which to test rho:
% (for convenience, just set s to be grid size)
s = 16;
rho1 = linspace(.1,1,s);
rho2 = linspace(.1,.9,s);

omega = [.3600 .5600];

res = grid_optimize_rho(dat,omega(1),omega(2),rho1,rho2);

min_mspe   = min(res.rho_mspe(:));
[row,col]  = find(res.rho_mspe == min_mspe);
best_rho = [res.rho1(row,col) res.rho2(row,col)];

res.rho_mspe
best_rho
min_mspe



%% Optimize omega
% Determine grid on which to test rho:
% (for convenience, just set s to be grid size)
s = 20;
omega1 = linspace(0.1,0.9,s);
omega2 = linspace(0.1,0.9,s);

rho = [.2600 .9550];

res = grid_optimize(dat,omega1,omega2,rho(1),rho(2));

min_mspe   = min(res.omega_mspe(:));
[row,col]  = find(res.omega_mspe == min_mspe);
best_omega = [res.omega1(row,col) res.omega2(row,col)];

res.omega_mspe
best_omega
min_mspe

% Current mspe 5.5896e-06


%% Optimize both
% Determine grid on which to test rho:
% (for convenience, just set s to be grid size)
s = 3;
omega1 = linspace(0.35,0.37,s);
omega2 = linspace(0.55,0.57,s);
rho1   = linspace(0.26,0.3,s);
rho2   = linspace(0.95,.955,s);
lambda = 1;

res = grid_optimize(dat,omega1,omega2,rho1,rho2);

min_mspe   = min(res.params_mspe(:));
[a,b,c,d]  = ind2sub(size(res.params_mspe),find(res.params_mspe == min_mspe));
best_omega = [res.omega1(a,b,c,d) res.omega2(a,b,c,d)];
best_rho = [res.rho1(a,b,c,d) res.rho2(a,b,c,d)];

res.params_mspe
best_omega
best_rho


% Current bests:
omega = [0.2 0.82];
rho   = [0.39 0.99];
lambda = 0.005;
% mspe = 5.026e-06

%% Optimize omega, rho and lambda
% Determine grid on which to test rho:
% (for convenience, just set s to be grid size)
s = 2;
omega1 = linspace(0.2,.22,s);
omega2 = linspace(0.82,.88,s);
rho1   = linspace(0.39,.42,s);
rho2   = linspace(0.985,.99,s);
lambda = linspace(0.001,.0015,s);

res = grid_optimize_all(dat,omega1,omega2,rho1,rho2,lambda);

min_mspe     = res.mspe;
%[a,b,c,d,e]  = ind2sub(size(res.params_mspe),find(res.params_mspe == min_mspe));
best_omega   = res.omega;
best_rho     = res.rho;
best_lambda  = res.lambda;

%res.params_mspe
best_omega
best_rho
best_lambda
min_mspe

% Outcome of optimization of each omega, rho on linspace(0.1,1,6)
% and lambda inspace (0.01,4,6): rho = omega = [0.82 0.82]
% lambda = 0.01, mspe = 5.6309e-06