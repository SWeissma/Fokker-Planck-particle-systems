%%% Toy example nonlinear 2 dimensional: Unimodal

clear;

%%% Setup

% nonlinear forward map
h = @(x) [x(2,:)*0.25+exp(-x(1,:))*(-0.25^2/2+0.25/2);x(2,:)*0.75+exp(-x(1,:))*(-0.75^2/2+0.75/2)];

% Parameter

N_x = 2;  % parameter space
R = 0.01*eye(N_x,N_x);  % noise covariance
m0 = [0;0]; % prior mean
P0 = 1*eye(N_x,N_x);  % prior covariance
K = 2;  % number of observations

% reference solution
% utrue = mvnrnd(m0,Gamma0,1)';   % drawn from prior
utrue = [0.0865; -0.8157];  % setup of the paper

% y = G(utrue)+chol(Gamma)'*randn(K,1);  % observed and perturbed by noise
y = [-0.0173;-0.573];   % setup of the paper


%%% Contourlines of the PDF

% prior
pi_0 = @(u) mvnpdf(u,m0,P0);

% potential
potential = @(u)1/2*norm((R^(1/2))\(y-h(u)))^2;

% posterior
pi_y = @(u)exp(-potential(u))*pi_0(u);

z_axis = linspace(-2,2,100)';
Cont = zeros(100,100);
for i1 = 1:100
    for i2 = 1:100
        Cont(i1,i2) = pi_y([z_axis(i1);z_axis(i2)]);
    end
end

figure(1)
clf(1)
contourf(z_axis,z_axis,Cont','LineWidth',2,'DisplayName','target distribution');hold on

save('reference_data.mat')

