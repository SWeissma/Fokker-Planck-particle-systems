%%% MH-MCMC 

clear;
load('reference_data.mat')

%%%%% MH-MCMC: gaussian-Kernel

% prior distribution: gaussian density
mu_0 = @(xi) mvnpdf(xi,zeros(I,1),Gamma0);

% potential
potential = @(xi)1/2*norm((Gamma^(1/2))\(y-G(xi)))^2;

% posterior distribution
mu_y = @(xi)exp(-potential(xi))*mu_0(xi);

% % Simulate MC until Q and drop 10*burn
q = 5;
burn = 10^(q-1);
Q = 2*10^q+10*burn;
U_MCMC = zeros(I,Q);

% Stepsize of the pCN-Kernel s
s = 0.07;

% Draw u_1 from mu_0 as initial value 
U_MCMC(:,1)=mvnrnd(m0,Gamma0,1);

% define acceptance probability
alpha = @(z1,z2) min(1,mu_y(z2)/mu_y(z1));

% acceptance rate (counter)
acc = 0;
tic;
for q = 1:Q-1
    % pCN proposal
    r = mvnrnd(sqrt(1-s^2)*U_MCMC(:,q),s^2*Gamma0)';
    v = rand;
    if v<= alpha(U_MCMC(:,q),r)
        U_MCMC(:,q+1)=r;
        acc = acc+1;
    else
        U_MCMC(:,q+1)=U_MCMC(:,q);
    end
    
end

runningtime_MCMC = toc;

% acceptance rate
acc = acc/(Q-1);

Xi_MCMC = U_MCMC(:,burn+1:end);
Uest5 = u(Xi_MCMC);
Uest_MCMC = mean(Uest5,2);

SD_MCMC = sqrt(var(Uest5')');
LB_MCMC = mean(Uest5,2)-SD_MCMC;
UB_MCMC = mean(Uest5,2)+SD_MCMC;

% save the results (might result in large file)
save('MCMC.mat','Uest_MCMC','SD_MCMC','LB_MCMC','UB_MCMC','acc','Uest5','runningtime_MCMC','Xi_MCMC','s');

fig4 = figure(4);    
clf(fig4)
set(fig4, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.3, 0.3]);

p3 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.7],'LineWidth',3,'DisplayName','MCMC');hold on
plot(xx_b(2:end-1),LB_MCMC,'-','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'-','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on
fill([xx',xx']',[LB_MCMC,UB_MCMC]','-','edgecolor',[0 0 1],'LineWidth',0.5);
p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on

title('(a) deterministic','FontSize',24)
xlabel('x')
legend('show',[p1,p3],'FontSize',16,'Location','southeast')
