% Simulations: Numerical integration of the 5D governing equations 
% of the fish dynamics
clear all, close all, clc  

%global variables linked to the function ODEfive
global alpha k rho  psi kappap kappav kappas

%non dimensional parameters
k=1; % lateral line gain
rho=0.1; %dipole length
alpha=0.18; %flow speed/fish speed ratio
psi=1; % self propulsion ratio v02/v01
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

% Initial conditions: 

%% TANDEM
lambda0 =0.8; %inter-fish distance (stream-wise)

ksi10=0+0.01*randn(1); %cross-stream coordinates
ksi20=0+0.01*randn(1);

thetaf10=pi+0.02*randn(1); %heading angles
thetaf20=pi+0.02*randn(1);

%% staggered 
% lambda0 =0.3; %inter-fish distance (stream-wise)
% 
% ksi10=0.15+0.01*randn(1); %cross-stream coordinates
% ksi20=-0.15+0.01*randn(1);
% 
% thetaf10=pi+0.02*randn(1); %heading anglese
% thetaf20=pi+0.02*randn(1);


%% ODE INTEGRATION using ode113 (better for non-linear systems)

y0 = [ksi10 ksi20 lambda0 thetaf10 thetaf20]'; %(5x1 vector with intial contidions)
tf=100; % time of integration (can be changed)
tspan = [0 tf]; 
opts=odeset('RelTol',1e-2,'AbsTol',1e-4,'Stats','on','OutputFcn',@odeplot);%set tolerances
[t,y] = ode113(@ODEfive,tspan,y0); %integration of the equations contained in ODEfive.m function

%% plot
figure(1) %cross-stream 
set(1,'units','normalized','name','variation on time ','numbertitle','off');
%figure('DefaultAxesFontSize',12)
plot(t,y(:,1),'r',t,y(:,2),'b','LineWidth',1.5)
xlabel('Time (s)','FontSize', 20),ylabel('\xi','FontSize', 20)
%title('\xi')
legend('\xi_1','\xi_2','location','best')

figure(2) % heading angles
set(2,'units','normalized','name','variation on time ','numbertitle','off');
plot(t,y(:,4),'g',t,y(:,5),'m','LineWidth',1.5)
%axis([0 tf 0 pi])
xlabel('Time (s)','FontSize', 20),ylabel('\theta','FontSize', 20)
legend('\theta_1','\theta_2','location','best')

figure(3) %stream-wise
set(3,'units','normalized','name','variation on time ','numbertitle','off');
plot(t,y(:,3),'b','LineWidth',1.5)
%axis([0 tf 0 pi])
xlabel('Time (s)','FontSize', 20),ylabel('\Lambda','FontSize', 20)
legend('\Lambda','location','best')

figure(4) %PHASE SPACE
set(4,'units','normalized','name','variation on time ','numbertitle','off');
plot(y(:,3),(y(:,1)-y(:,2)),'m','LineWidth',1.5)
%axis([0 tf 0 pi])
xlabel('\Lambda','FontSize', 20),ylabel('\Delta \xi','FontSize', 20)
