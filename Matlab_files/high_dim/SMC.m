%%% SMC combined with Langevin dynamics

clear;
load('reference_data.mat');

N = 250; % amount of sequential steps

% prior distribution: gaussian density
pi_0 = @(xi) exp(-1/2*sum((sqrtm(Gamma0)\(xi-m0)).^2,1));

% potential tempered
likelihood_N = @(xi)exp(-1/2*1/N*sum((sqrtm(Gamma)\(G(xi)-y)).^2,1));

% initialization
M = 512;
Xi_SMC = mvnrnd(m0,Gamma0,M)';  % initial ensemble

% initial weights for SMC
weights = 1/M*ones(1,M);    % initial weights
M_tol = 100;    % tolerance value for ESS

% acceptance probability
alpha = @(z1,z2) min(1,(likelihood_N(z2).*pi_0(z2))./(likelihood_N(z1).*pi_0(z1)));

% initialize time and spread of the particle system
time_SMC = 0;
spread_SMC = 0;

tic;
% run the sequential Langevin method
for k = 1:N
    
    % update weights and compute ESS
    weights = weights.*likelihood_N(Xi_SMC);
    weights = weights/sum(weights);
    ESS = sum(weights)^2/sum(weights.^2);
    
    % running time for Langevin dynamics
    End_time = 0.004;   
    
    % run Langevin dynamic
    [Xi_path2,time2] = Langevin_SMC(G,Gamma0,Xi_SMC,y,End_time,Gamma,k/N);
    Xi_SMC = reshape(Xi_path2(end,:),[I,M]);
    
    % collect time
    time_SMC = [time_SMC, time_SMC(end)+time2]; 
    
    % compute the spread of the particles
    spread2 = zeros(1,length(time2));
    for tk = 1:length(time2)
        Xihelp = reshape(Xi_path2(tk,:),[I,M]);
        spread2(tk) = 1/M*norm(Xihelp-mean(Xihelp,2),'fro')^2;
    end
    spread_SMC = [spread_SMC,spread2];
    
    % if ESS <= M_tol: resample according to weights and reset weights
    if ESS <= M_tol  
        resampling = rand(M,1);
        Q = cumsum(weights).*ones(M,M)-resampling;
        Ind = (Q<0);
        i = sum(Ind,2)+1;
        Xi_SMC = Xi_SMC(:,i);
        weights = 1/M*ones(1,M);
    end
    
    % update acceptance probability
    alpha = @(z1,z2) min(1,(likelihood_N(z2).^(k+1).*pi_0(z2))./(likelihood_N(z1).^(k+1).*pi_0(z1)));
    
end
time_SMC(1) =[];
spread_SMC(1) =[];

% measure final running time
runningtime_SMC = toc;


% Parameter estimate
Uest_SMC = mean(u(Xi_SMC),2);
SD_SMC = sqrt(var(u(Xi_SMC)')');
LB_SMC = mean(u(Xi_SMC),2)-SD_SMC;
UB_SMC = mean(u(Xi_SMC),2)+SD_SMC;
   

% load MCMC results to compare
load('MCMC.mat')


%%% Parameter estimation

fig3 = figure(3);    
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.3, 0.3]);

plot(xx_b(2:end-1),LB_SMC,':','color',[1 0.5 0.3],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_SMC,':','color',[1 0.5 0.3],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),LB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
fill([xx',xx']',[LB_SMC,UB_SMC]',':','edgecolor',[0.8 0.3 0.1],'LineWidth',0.5);
fill([xx',xx']',[LB_MCMC,UB_MCMC]','--','edgecolor',[0 0 1],'LineWidth',0.5);

p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on
p3 = plot(xx_b,[0;Uest_SMC;0],'--','Color',[0.8 0.3 0.1],'LineWidth',3,'DisplayName','SMC');hold on
p4 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.5],'LineWidth',3,'DisplayName','MCMC');hold on

title('(b) SMC','FontSize',24)
xlabel('x')
legend('show',[p1,p3,p4],'FontSize',16,'Location','northwest')

%%% scatter plots

fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.6]);

subplot(2,2,3)
s1 = scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
s2 = scatter(Xi_SMC(1,:),Xi_SMC(2,:),'MarkerFaceColor',[0.8 0.3 0.1],'MarkerEdgeColor',[0.3 0 0],'DisplayName','SMC');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(c) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

subplot(2,2,4)
s4 = scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
s5 = scatter(Xi_SMC(1,:),Xi_SMC(3,:),'MarkerFaceColor',[0.8 0.3 0.1],'MarkerEdgeColor',[0.3 0 0],'DisplayName','SMC');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(d) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);



% nodes of time for running Langevin 
TT = time_SMC;

% Parameters
M2 = M;   % ensemble size   

% prior sample
u0 = mvnrnd(m0,Gamma0,M2)';  % initial ensemble

tic;
% Final time for Langevin dynamics
End_time = N*End_time;  
% run Langevin dynamics
[Xi_path2,time2] = Langevin_correction_compare(G,Gamma0,u0,y,End_time,Gamma,TT);

% measure final running time
runningtime2 = toc;

% compute the spread of the particles
spread2 = zeros(1,length(time2));
for tk = 1:length(time2)
    Xihelp = reshape(Xi_path2(tk,:),[I,M2]);
    spread2(tk) = 1/M2*norm(Xihelp-mean(Xihelp,2),'fro')^2;
end
Xi_sample2 = reshape(Xi_path2(end,:),[I,M2]);


%%% scatter plots
subplot(2,2,1)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(2,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','Langevin');hold on

xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(a) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

subplot(2,2,2)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(3,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','Langevin');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(b) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)



Uest5 = u(Xi_MCMC(:,end-20*9999:20:end));
Uest_MCMC = mean(Uest5,2);

SD_MCMC = sqrt(var(Uest5')');
LB_MCMC = mean(Uest5,2)-SD_MCMC;
UB_MCMC = mean(Uest5,2)+SD_MCMC;


%%% parameter estimate
fig6 = figure(6);   
clf(fig6)
set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);
subplot(1,2,1)

Uest_sample2 = u(Xi_sample2);
Uest_sample_mean2 = mean(Uest_sample2,2);
SD_sample2 = sqrt(var(Uest_sample2')');
LB_sample2 = Uest_sample_mean2-SD_sample2;
UB_sample2 = Uest_sample_mean2+SD_sample2;

plot(xx,LB_sample2,':','color',[0 1 0],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx,UB_sample2,':','color',[0 1 0],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),LB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
fill([xx;xx],[LB_sample2';UB_sample2'],':','edgecolor',[0 1 0],'LineWidth',0.5);
fill([xx',xx']',[LB_MCMC,UB_MCMC]','--','edgecolor',[0 0 1],'LineWidth',0.5)


p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on
p2 = plot(xx_b,[0;Uest_sample_mean2;0],':','Color',[0 0.5 0],'LineWidth',3,'DisplayName','estimate');hold on
p4 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.5],'LineWidth',3,'DisplayName','MCMC');hold on

title('(a) stochastic','FontSize',24)
xlabel('x')
legend('show',[p1,p2,p4],'FontSize',16,'Location','northwest')

subplot(1,2,2)
plot(xx_b(2:end-1),LB_SMC,':','color',[1 0.5 0.3],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_SMC,':','color',[1 0.5 0.3],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),LB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
fill([xx',xx']',[LB_SMC,UB_SMC]',':','edgecolor',[0.8 0.3 0.1],'LineWidth',0.5);
fill([xx',xx']',[LB_MCMC,UB_MCMC]','--','edgecolor',[0 0 1],'LineWidth',0.5);


p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on
p3 = plot(xx_b,[0;Uest_SMC;0],':','Color',[0.8 0.3 0.1],'LineWidth',3,'DisplayName','estimate');hold on
p4 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.5],'LineWidth',3,'DisplayName','MCMC');hold on

title('(b) SMC','FontSize',24)
xlabel('x')
legend('show',[p1,p3,p4],'FontSize',16,'Location','northwest')


%%% compare spread of particle systems
fig7 = figure(7);
clf(fig7)
set(fig7, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);

plot(time2,spread2,'-','Color',[0 0 0.7],'LineWidth',2,'DisplayName','Langevin');hold on
plot(time_SMC,spread_SMC,':','Color',[0.7 0 0],'LineWidth',2,'DisplayName','SMC');hold on
ylim([0, 0.5])
xlabel('time','FontSize',20)
l = legend('show','Location','northeast');
set(l,'Interpreter','latex','FontSize',20);
title('(a) ensemble spread','FontSize',24)


