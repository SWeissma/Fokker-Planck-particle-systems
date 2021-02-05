%%% Interacting Langevin particle dynamics

clear;
% load the data
load('reference_data.mat')

% path to functions
path('./Functions',path);

% Parameters
M2 = 200;   % ensemble size   

gam = 1;    % scaling parameter of localization

% prior sample
x0 = mvnrnd(m0,P0,M2)';  % initial ensemble
x0 = x0(:);


% define drift and diffusion
drift = @(t,U) Langevin_drift(U,N_x,h,R,P0,m0,y);
diff = @(t,U) Langevin_diff(U,N_x);
obj = sde(drift,diff,'StartState',x0);

tic;
nPeriods = 10000;   % number of iterations
dt = 1/1000;    % stepsize
[Xest2_path,time2] = simByEuler(obj,nPeriods,'DeltaTime',dt);

runningtime = toc;


save('SDE_Langevin.mat')
