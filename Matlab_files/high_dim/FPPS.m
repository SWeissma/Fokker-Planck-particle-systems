%%% Fokker--Planck Particle System

clear;

% load the data
load('reference_data.mat')

% path to functions
path('./Functions',path);

% Parameters
M = 2^7;   % ensemble size 
% M = 2^8;
% M = 2^9;

% prior sample
xi0 = mvnrnd(m0,Gamma0,M)';  % initial ensemble
xi0_quer = mean(xi0,2);   % initial ensemble mean

% product kernel
delta = 1/(I+4);

tic;

% options for the ode setup
tspan = [0,1];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% run until T=1;
[t, Xi] = ode45(@(t,Xi) odesystem_FPPS(Xi,G,y,Gamma,I,Gamma0,m0,delta), tspan, xi0(:),options);
xi1 = Xi(end,:);
Xi1 = reshape(Xi(end,:),[I,M]);
xi1_quer = mean(Xi1,2);
Pxx = (1/M)*Xi1*Xi1'-xi1_quer*xi1_quer';
B_est = diag(Pxx);
B = (4/((2+I)*M))^delta*diag(B_est);

% run until T = 100 with fixed kernel
tspan = [0,100];
[t, Xi] = ode45(@(t,Xi) odesystem_FPPS_fix(Xi,G,y,Gamma,I,Gamma0,m0,B), tspan, xi0(:),options);

runningtime1 = toc;

Xi_est = reshape(Xi(end,:),[I,M]);
xi_quer = mean(Xi_est,2);


% sampling and resampling of the kernel density estimate
Z = sampling_kernel(Xi_est,B,100);
Xi_sample = resampling(Z,Xi_est,B,512);
Uest_sample = u(Xi_sample);

% save the results
save('FPPS_M128.mat','Xi_sample','Uest_sample','Xi','t','G','y','Gamma','I','Gamma0','m0','B','delta','xi0','M','runningtime1')
% save('FPPS_M256.mat','Xi_sample','Uest_sample','Xi','t','G','y','Gamma','I','Gamma0','m0','B','delta','xi0','M','runningtime1')
% save('FPPS_M512.mat','Xi_sample','Uest_sample','Xi','t','G','y','Gamma','I','Gamma0','m0','B','delta','xi0','M','runningtime1')





