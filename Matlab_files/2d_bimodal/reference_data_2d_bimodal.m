%%% Toy example nonlinear 2 dimensional: Bimodal

clear;
close all;

%%% Setup

% nonlinear forward map
h = @(u) (u(1,:)-u(2,:)).^2;

% gradient of the forward map
gradh = @(u) [2*(u(1,:)-u(2,:));-2*(u(1,:)-u(2,:))];

% Parameter

N_x = 2;  % parameter space
R = 1;    % noise covariance
m0 = [0;0];    % prior mean
P0 = eye(N_x,N_x);  % prior covariance

% reference solution
% utrue = mvnrnd(m0,Gamma0,1)';   % drawn from prior
xtrue = [-1.5621; -0.0021];     % setup of the paper

% y = G(utrue)+chol(Gamma)'*randn(1);  % observed and perturbed by noise
y = 4.2297;  % setup of the paper


%%% Contourlines of the PDF

% prior
pi_0 = @(u) mvnpdf(u,m0,P0);
% potential
potential = @(u)1/2*norm((R^(1/2))\(y-h(u)))^2;
% posterior
pi_y = @(u)exp(-potential(u))*pi_0(u);

% Contourlines
z_axis = linspace(-4,4,100)';
Cont = zeros(100,100);
for i1 = 1:100
    for i2 = 1:100
        Cont(i1,i2) = pi_y([z_axis(i1);z_axis(i2)]);
    end
end

figure(1)
clf(1)
contourf(z_axis,z_axis,Cont','LineWidth',2,'DisplayName','target distribution');
legend('show','Location','southwest','FontSize',16)

% save the data
save('reference_data.mat')

