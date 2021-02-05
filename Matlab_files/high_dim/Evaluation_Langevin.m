
clear;
close all;

% load the corresponding data:
% note that one has to run Langevin.m and 
% Langevin_corrected.m for M2 = 16,32,64 first.
load('reference_data.mat');
load('MCMC.mat')

%%%%%%%%%%%%%%%%%%% Evaluation Langevin without correction
%%% scatter plots

load('Langevin_without_correction_M16.mat');

fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.6]);

subplot(2,2,1)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(2,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','M = 16');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(a) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

subplot(2,2,2)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(3,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','M = 16');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(b) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

load('Langevin_without_correction_M64.mat');

subplot(2,2,3)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(2,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','M = 64');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(c) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

subplot(2,2,4)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(3,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','M = 64');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(d) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)


%%% Plotting spread

fig3 = figure(3);
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);
subplot(1,2,1)
load('Langevin_without_correction_M16.mat');
plot(time(1:end-1),spread(1:end-1),'-','Color',[0 0 0.4],'LineWidth',2,'DisplayName','M = 16');hold on
load('Langevin_without_correction_M32.mat');
plot(time(1:end-1),spread(1:end-1),'-','Color',[0 0 0.6],'LineWidth',2,'DisplayName','M = 32');hold on
load('Langevin_without_correction_M64.mat');
plot(time(1:end-1),spread(1:end-1),'-','Color',[0 0 0.8],'LineWidth',2,'DisplayName','M = 64');hold on
xlabel('time','FontSize',20)
ylim([0, 0.5])
l = legend('show','Location','northeast');
set(l,'Interpreter','latex','FontSize',20);
title('(a) ensemble spread','FontSize',24)


%%%%%%%%%%%%%%%%%%% Evaluation Langevin with correction
%%% scatter plots

load('Langevin_M16.mat');

fig4 = figure(4);
clf(fig4)
set(fig4, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.6]);

subplot(2,2,1)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(2,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','corrected, M = 16');hold on

xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(a) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

subplot(2,2,2)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(3,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','corrected, M = 16');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(b) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

load('Langevin_M64.mat');

subplot(2,2,3)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(2,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','corrected, M = 64');hold on

xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(c) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)

subplot(2,2,4)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample2(1,:),Xi_sample2(3,:),'MarkerFaceColor',[0 .6 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','corrected, M = 64');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(d) particle system');
title(str,'FontSize',20)
legend('show','Location','southwest','FontSize',16)


%%% Plotting spread and residuals

figure(3);
subplot(1,2,2)
load('Langevin_M16.mat');
plot(time2(1:end-1),spread2(1:end-1),'-','Color',[0 0 0.4],'LineWidth',2,'DisplayName','M = 16');hold on
load('Langevin_M32.mat');
plot(time2(1:end-1),spread2(1:end-1),'-','Color',[0 0 0.6],'LineWidth',2,'DisplayName','M = 32');hold on
load('Langevin_M64.mat');
plot(time2(1:end-1),spread2(1:end-1),'-','Color',[0 0 0.8],'LineWidth',2,'DisplayName','M = 64');hold on
ylim([0, 0.5])
xlabel('time','FontSize',20)
l = legend('show','Location','northeast');
set(l,'Interpreter','latex','FontSize',20);
title('(b) ensemble spread (corrected)','FontSize',24)



%%% Parameter estimation


fig6 = figure(6);   
load('Langevin_M64.mat');

Uest5 = u(Xi_MCMC(:,end-20*9999:20:end));
Uest_MCMC = mean(Uest5,2);

SD_MCMC = sqrt(var(Uest5')');
LB_MCMC = mean(Uest5,2)-SD_MCMC;
UB_MCMC = mean(Uest5,2)+SD_MCMC;

clf(fig6)
set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);
subplot(1,2,2)

p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on
p3 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.7],'LineWidth',3,'DisplayName','MCMC');hold on
plot(xx_b(2:end-1),LB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on

Uest_sample2 = u(Xi_sample2);
Uest_sample_mean2 = mean(Uest_sample2,2);
SD_sample2 = sqrt(var(Uest_sample2')');
LB_sample2 = Uest_sample_mean2-SD_sample2;
UB_sample2 = Uest_sample_mean2+SD_sample2;
p2 = plot(xx_b,[0;Uest_sample_mean2;0],':','Color',[0 0.7 0],'LineWidth',3,'DisplayName','derivative-free');hold on
plot(xx,LB_sample2,':','color',[0 1 0],'LineWidth',3,'DisplayName','credible set');hold on
plot(xx,UB_sample2,':','color',[0 1 0],'LineWidth',3,'DisplayName','credible set');hold on
fill([xx;xx],[LB_sample2';UB_sample2'],':','edgecolor',[0 1 0],'LineWidth',0.5);
fill([xx',xx']',[LB_MCMC,UB_MCMC]','--','edgecolor',[0 0 1],'LineWidth',0.5);



title('(b) stochastic','FontSize',24)
xlabel('x')
legend('show',[p1,p2,p3],'FontSize',16,'Location','northwest')

