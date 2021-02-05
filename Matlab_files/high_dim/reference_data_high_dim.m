%%% Toy example nonlinear high dimensional:

clear;
close all;

% Parameter
K = 2^5-1;      % number of observations
l = 8;          % discretization of the PDE
I_u = 2^l;  
xx_b = linspace(0,1,I_u+1);     % Grid of the PDE with boundary
xx = xx_b(2:end-1);     % Grid of the PDE without boundary

% KL Basis
V = zeros(I_u-1,I_u-1);
D = zeros(I_u-1,I_u-1);
for i = 1:(I_u-1)
    V(:,i) = sqrt(2)*sin(i*pi*xx)';
    D(i,i) = (i^2)^(-1.5);
end

% truncating KL expansion by I
I = 2^5;    

% evaluation of the KL expansion
u = @(xi) V(:,1:I)*xi;

% Observation operator
Obs_mat = observation_matrix(K,I_u);

% forward map
G = @(xi) Obs_mat*pdesolvenonl(l,u(xi));
G_u =@(uu) Obs_mat*pdesolvenonl(l,uu);

% parameter prior and posterior
m0 = zeros(I,1);                % prior mean
Gamma0 = eye(I,I)*D(1:I,1:I);     % prior covariance

% compute a reference solution by truncating KL expansion by 2^2
Xi_true = sqrt(D(1:I,1:I))*randn(I,1);
Xi_true(2^2+1:end) = zeros(I-2^2,1);

% evaluate KL expansion
utrue = V(:,1:I)*Xi_true;
% noise covariance
Gamma = 0.0001*eye(K,K);        

% compute noisy observation
y = G_u(utrue) + 0.1*sqrt(Gamma)*randn(K,1);

% plot the parameter and the corresponding observation
figure(3)
clf(3)
p1 = plot(xx_b,[0;utrue;0],'black','LineWidth',2,'DisplayName','truth');
figure(4)
clf(4)
p2 = plot(linspace(0,1,length(y)+2),[0;y;0],'+');

% save the data
save('reference_data.mat')