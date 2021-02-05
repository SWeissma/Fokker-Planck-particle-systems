%%% Fokker--Planck Particle System

clear;

% load the data
load('reference_data.mat')

% path to functions
path('./Functions',path);

% Parameters
M = 200;   % ensemble size   

% prior sample
x0 = mvnrnd(m0,P0,M)';  % initial ensemble
x0_quer = mean(x0,2);   % initial ensemble mean

% construction of gaussian kernel
Pxx0 = 1/M*x0*x0'-x0_quer*x0_quer';   % sample covariance
alpha = 0.05;    % kernel scaling
B = alpha*Pxx0;

tic;
% options for the ode setup
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% FPPS derivative-free
[t, X] = ode45(@(t,U) odesystem_FPPS(U,h,y,R,N_x,P0,m0,B),[0 1],x0,options);
Xest = reshape(X(end,:),[N_x,M]);

runningtime = toc;

save('ode_FPPS.mat','X','t','h','y','R','N_x','P0','m0','B','alpha','x0','M','runningtime')
