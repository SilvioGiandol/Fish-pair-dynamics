%Heatmaps generations TANDEM, stable and unstable solutions function of
%lambda and alpha
clear all, clc, close all
global r kappa a l  ps kappas kappap kappav

ranger= 30; % number of stepsize for linespace (best if 30)

%% K=0 RHO=0.1
%non dimensional parameters
kappa=0; %lateral line
r=0.1; %dipole length rho
BL=r/0.2;
a=0.18; %flow speed/fish speed ratio
ps=1; % self propulsion ratio
kappap= 0; %0.03 attraction gain
kappav=0; %0.558 alignment gain
kappas=0; %alignment gain due to surroundings

in1=0; %for tandem
X1=[in1 pi pi]; % initial conditions for numerical solver fsolve [ksi theta1 theta2] 
omega1=zeros(length(a),length(l)); % matrix of natural frequencies
omega2=zeros(length(a),length(l));
eigenva= zeros(length(a),length(l)); %size (a) RESET
g=1; %adder columns
for a=linspace(0,0.2,ranger) %range of flow velocity
    f=1; %adder rows
    for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L)
        lam=l;
        opts.Algorithm = 'levenberg-marquardt';
        Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
        zita1 = Root1(1); 
        th1=Root1(2);
        th2=Root1(3);
        x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and check equations
        E=equationmap(x) %function to check equations are satisfied
        [EIGENMAP,om1,om2]=jacob(x); % check jacobian 
        if E==1 && abs(zita1)<0.0001 %condition for tandem
            if EIGENMAP==1  %negative (STABLE SPIRAL)
                eigenva(f,g)=-1;
                omega1(f,g)= NaN;
                omega2(f,g) = NaN;
            elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
                eigenva(f,g)=0;
                omega1(f,g)= om1;
                omega2(f,g) = om2;
            else %positive  (UNSTABLE)
                eigenva(f,g)=1;
                omega1(f,g)= NaN;
                omega2(f,g) = NaN;
            end
        else
            eigenva(f,g)=1;
            omega1(f,g)= NaN;
            omega2(f,g) = NaN;
        end
        f=f+1;
    end
    g=g+1;
end

%%HEATMAP GENERATION
a= linspace(0,0.2,ranger); % range of velocity (x axis)
l= linspace((10^-16)/BL,1/BL,ranger);%TANDEM %range of inter-fish (y-axis) 
[velocity,interfish]= meshgrid(a,l); % mesh

%colors:
cMap = [0, 0, 1; ...   % Blue for -1 positive 
  0, 1, 1; ...       % Greenfor 0 zero
  1, 0, 0]     % Red for 1 negative

F1=figure(1);

[c,h]=contourf(velocity, interfish,eigenva,50);
set(h, 'edgecolor','none');
xlim('auto')
ylim('auto')
colormap(cMap);
cb=colorbar;    
ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
caxis([-1 1]);  % axis of the colorbar
xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
hold on;
h=3.1/r  %mm because 3.1 is the dipole length in mm

%IF we want to save the graphs:
% saveas(F1, 'tandemK0R01.png')
% close(F1) %close the graph
%if we want to save the frequencies:
% save('tandemomega1k0r01.mat','omega1'); %save omega1
% save('tandemomega2k0r01.mat','omega2');  %save omega2


%% ALL THE OTHER CASES 
% %% K=0 RHO=0.05
% r=0.05;
% in1=0; %for tandem
% X1=[in1 pi pi]; 
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     a
%     f=1; %adder rows
%     for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && abs(zita1)<0.0001 %condition for tandem
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% 
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,1/0.5,ranger);%TANDEM
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% F2=figure(2);
% [c,h]=contourf(velocity, interfish,eigenva,50);
% set(h, 'edgecolor','none');
% xlim('auto')
% ylim('auto')
% colormap(cMap);
% cb=colorbar;    
% ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
% caxis([-1 1]);  % axis of the colorbar
% xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
% hold on;
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F2, 'tandemK0R005.png')
% close(F2)
% %% K=0 RHO=0
% r=0; %infinity channel
% in1=0; %for tandem
% X1=[in1 pi pi]; 
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     a
%     f=1; %adder rows
%     for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=eigenmap(x); % check jacobian with the values 
%         if E==1 && abs(zita1)<0.0001 %condition for tandem
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% 
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,1/0.5,ranger);%TANDEM
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% F3=figure(3);
% [c,h]=contourf(velocity, interfish,eigenva,50);
% set(h, 'edgecolor','none');
% xlim('auto')
% ylim('auto')
% colormap(cMap);
% cb=colorbar;    
% ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
% caxis([-1 1]);  % axis of the colorbar
% xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
% hold on;
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F3, 'tandemK0R0.png')
% close(F3)
% %% K=1 RHO=0.1
% kappa=1;  % LATERAL-LINE SENSING
% r=0.1;
% kappap=0; %0.03 attraction gain
% kappav=0; %0.558 allignment gain
% kappas=0;
% ps=1;
% in1=0; %for tandem
% X1=[in1 pi pi]; 
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     a
%     f=1; %adder rows
%     for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && abs(zita1)<0.0001 %condition for tandem
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% 
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,1/0.5,ranger);%TANDEM
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% F4=figure(4);
% [c,h]=contourf(velocity, interfish,eigenva,50);
% set(h, 'edgecolor','none');
% xlim('auto')
% ylim('auto')
% colormap(cMap);
% cb=colorbar;    
% ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
% caxis([-1 1]);  % axis of the colorbar
% xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
% hold on;
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F4, 'tandemK1R01.png')
% close(F4)
% %% K=1 RHO=0.05
% r=0.05;
% in1=0; %for tandem
% X1=[in1 pi pi]; 
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     a
%     f=1; %adder rows
%     for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && abs(zita1)<0.0001 %condition for tandem
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% 
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,1/0.5,ranger);%TANDEM
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% F5=figure(5);
% [c,h]=contourf(velocity, interfish,eigenva,50);
% set(h, 'edgecolor','none');
% xlim('auto')
% ylim('auto')
% colormap(cMap);
% cb=colorbar;    
% ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
% caxis([-1 1]);  % axis of the colorbar
% xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
% hold on;
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F5, 'tandemK1R005.png')
% close(F5)
% %% K=0 RHO=0
% r=0; %infinity channel
% in1=0; %for tandem
% X1=[in1 pi pi]; 
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     a
%     f=1; %adder rows
%     for l=linspace(10^-16,1,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X1,opts) %TANDEM
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && abs(zita1)<0.0001 %condition for tandem
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% 
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,1/0.5,ranger);%TANDEM
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% F6=figure(6);
% [c,h]=contourf(velocity, interfish,eigenva,50);
% set(h, 'edgecolor','none');
% xlim('auto')
% ylim('auto')
% colormap(cMap);
% cb=colorbar;    
% ylabel(cb, 'Max real part of eigenvalues','FontSize', 14);
% caxis([-1 1]);  % axis of the colorbar
% xlabel('\alpha','FontSize', 14),ylabel('\Lambda (BL)','FontSize', 14)
% hold on;
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F6, 'tandemK1R0.png')
% close(F6)