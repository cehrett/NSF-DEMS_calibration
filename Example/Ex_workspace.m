% Simulation workspace

clc; clear all; close all;

%% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);

%% Take a look at surfaces of example function output
theta1=linspace(0,6);
theta2=linspace(0,6);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 T2 ] = meshgrid(theta1,theta2); 
c = cat(2,T1',T2');
Theta = reshape(c,[],2);

% Set other parameters
alpha = repmat(2,size(Theta,1),1);
c     = alpha;

% Now get output values
opc = Ex_sim(Theta,alpha,c);

% Now take a look at the surfaces
oscls = reshape(opc(:,1),[],length(theta1));
perfs = reshape(opc(:,2),[],length(theta1));
costs = reshape(opc(:,3),[],length(theta1));

surf(theta1,theta2,oscls,'FaceColor','red','EdgeColor','flat');
axis vis3d;

hold on;
mult = max(oscls)/max(perfs);
surf(theta1,theta2,mult*perfs,'FaceColor','blue','EdgeColor','flat'); 
axis vis3d;

mult = max(oscls)/max(costs);
surf(theta1,theta2,mult*costs,'FaceColor','green','EdgeColor','flat'); 
axis vis3d;
xlabel('\theta_2'); ylabel('\theta_1'); zlabel('Outcomes');

