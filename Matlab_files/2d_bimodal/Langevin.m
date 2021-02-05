%%% Interacting Langevin particle dynamics

clear;
% load the data
load('reference_data.mat')

% path to functions
path('./Functions',path);

% Parameters
M2 = 200;   % ensemble size   

gam = 0.5;    % scaling parameter of localization

% prior sample
x0 = mvnrnd(m0,P0,M2)';  % initial ensemble
x0 = x0(:);

% define drift and diffusion
drift = @(t,U) Langevin_drift(U,N_x,h,R,P0,m0,y);
diff = @(t,U) Langevin_diff(U,N_x);
obj = sde(drift,diff,'StartState',x0);

drift_loc = @(t,U) Langevin_drift_loc(U,N_x,h,R,P0,m0,y,gam);
diff_loc = @(t,U) Langevin_diff_loc(U,N_x,gam);
obj_loc = sde(drift_loc,diff_loc,'StartState',x0);

drift_exact = @(t,U) Langevin_drift_exact(U,N_x,h,R,P0,m0,y,gradh);
diff_exact = @(t,U) Langevin_diff_exact(U,N_x);
obj_exact = sde(drift_exact,diff_exact,'StartState',x0);


tic;
nPeriods = 1000;   % number of iterations
dt = 1/100;    % stepsize

% Langevin derivative-free
[Xest2_path,time2] = simByEuler(obj,nPeriods,'DeltaTime',dt);

% Langevin exact
[Xest2_exact_path,time2_exact] = simByEuler(obj_exact,nPeriods,'DeltaTime',dt);

% Langevin derivative-free localised
[Xest2_loc_path,time2_loc] = simByEuler(obj_loc,nPeriods,'DeltaTime',dt);

runningtime = toc;

% save the results
save('SDE_Langevin.mat')
