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
Pxx0 = (1/M)*x0*x0'-x0_quer*x0_quer';   % sample covariance
alpha = 0.01;    % kernel scaling
B = alpha*Pxx0;

% localisation scaling
gam = 0.5;

tic;
% options for the ode setup
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% FPPS derivative-free
[t, X] = ode45(@(t,U) odesystem_FPPS(U,h,y,R,N_x,P0,m0,B),[0 10],x0,options);
Xest = reshape(X(end,:),[N_x,M]);

% FPPS exact
[t2, X2] = ode45(@(t2,U2) odesystem_FPPS_exact(U2,h,y,R,N_x,P0,m0,B,gradh),[0 10],x0,options);
Xest_d = reshape(X2(end,:),[N_x,M]);

% FPPS derivative-free localised
[t3, X3] = ode45(@(t3,U3) odesystem_FPPS_loc(U3,h,y,R,N_x,P0,m0,B,gam),[0 10],x0,options);
Xest_loc = reshape(X3(end,:),[N_x,M]);

runningtime = toc;

% save the results
save('ode_FPPS.mat','X','X2','X3','t','t2','t3','h','y','R','N_x','P0','m0','B','gam','alpha','x0','M','runningtime')
