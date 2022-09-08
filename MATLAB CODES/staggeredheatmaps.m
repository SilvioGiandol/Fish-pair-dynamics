%Heatmaps generations STAGGERED
clear all, clc, close all
global r kappa a l  ps kappas kappap kappav
ranger= 30; % number of stepsize for linespace

%% IS THE SAME AS TANDEMHEATMAPS
%% K=0 RHO=0.1
kappa=1;
r=0.1;
in2=0.09; %for staggere T-dipole  %0.25 for A-dipole
X2=[in2 pi pi];  
kappap=0; %0.03 attraction gain
kappav=0; %0.558 allignment gain
kappas=0; %alignmnent gain
ps=1;
R2=zeros(length(a),length(l));
omega1=zeros(length(a),length(l));
omega2=zeros(length(a),length(l));
eigenva= zeros(length(a),length(l)); %size (a) RESET
g=1; %adder columns
for a=linspace(0,0.2,ranger) %range of velocity
    f=1; %adder rows
    for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
        lam=l;
        opts.Algorithm = 'levenberg-marquardt';
        Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
        zita1 = Root1(1); %crea il vettore di root 1
        th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
        th2=Root1(3);
        x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
        E=equationmap(x)
        [EIGENMAP,om1,om2]=jacob(x); % check jacobian with the values 
        if E==1 && zita1>0.0001  %condition for side-by-side
            if EIGENMAP==1  %negative (STABLE SPIRAL)
                eigenva(f,g)=-1;
                R2(f,g) = NaN;
                omega1(f,g)= NaN;
                omega2(f,g) = NaN;
            elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
                eigenva(f,g)=0;
                R2(f,g) = 2*(abs(Root1(1)));
                omega1(f,g)= om1;
                omega2(f,g) = om2;
            else %positive  (UNSTABLE)
                eigenva(f,g)=1;
                R2(f,g) = NaN;
                omega1(f,g)= NaN;
                omega2(f,g) = NaN;
            end
        else
            eigenva(f,g)=1;
            R2(f,g) = NaN;
            omega1(f,g)= NaN;
            omega2(f,g) = NaN;
        end
        f=f+1;
    end
    g=g+1;
end
a= linspace(0,0.2,ranger); % range of velocity (x axis)
l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
[velocity,interfish]= meshgrid(a,l);
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
%Show the height of the channel
h=3.1/r  %mm because 3.1 is the dipole length in mm
%saveas(F1, 'staggeredK0R01.png')
%close(F1)
%save('k0r01.mat','R2'); %save delta ksi
save('omega1k0r01.mat','omega1'); %save omega1
save('omega2k0r01.mat','omega2');  %save omega2

% %% K=0 RHO=0.05
% r=0.05;  
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% R2=zeros(length(a),length(l));
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     f=1; %adder rows
%     for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l;
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && zita1>0.0001  %condition for side-by-side
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%                 R2(f,g) = NaN;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%                 R2(f,g) = 2*(abs(Root1(1)));
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%                 R2(f,g) = NaN;
%             end
%         else
%             eigenva(f,g)=1; %THIS HAS TO BE NAN
%             R2(f,g) = NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% 
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
% Show the height of the channel
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F2, 'staggeredK0R005.png')
% close(F2)
% save('k0r005.mat','R2');
% %% K=0 RHO=0
% r=0;
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     f=1; %adder rows
%     for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l;
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && zita1>0.0001  %condition for side-by-side
%             R2(f,g) = Root1(1);
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%             R2(f,g) = NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% 
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
% %Show the height of the channel
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F3, 'staggeredK0R0.png')
% close(F3)
% % K=1 RHO=0.1
% kappa=1;
% r=0.1;
% kappap=0; %0.03 attraction gain
% kappav=0; %0.558 allignment gain
% kappas=0;
% ps=1;
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% R2=zeros(length(a),length(l));
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     f=1; %adder rows
%     for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l;
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && zita1>0.0001  %condition for side-by-side
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%                 R2(f,g) = NaN;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%                 R2(f,g) = 2*(abs(Root1(1)));
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%                 R2(f,g) = NaN;
%             end
%         else
%             eigenva(f,g)=1;
%             R2(f,g) = NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% 
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
% Show the height of the channel
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F4, 'staggeredK1R01.png')
% close(F4)
% save('k1r01.mat','R2');
% % K=1 RHO=0.05
% r=0.05;  
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% R2=zeros(length(a),length(l));
% g=1; %adder columns
% 
% for a=linspace(0,0.2,ranger) %range of velocity
%     f=1; %adder rows
%     for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l;
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && zita1>0.0001  %condition for side-by-side
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%                 R2(f,g) = NaN;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%                 R2(f,g) = 2*(abs(Root1(1)));
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%                 R2(f,g) = NaN;
%             end
%         else
%             eigenva(f,g)=1;
%             R2(f,g) = NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% 
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
% Show the height of the channel
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F5, 'staggeredK1R005.png')
% close(F5)
% save('k1r005.mat','R2');
% %% K=1 RHO=0
% r=0;
% in2=0.12; %for staggere T-dipole
% X2=[in2 pi pi];  
% eigenva= zeros(length(a),length(l)); %size (a) RESET
% g=1; %adder columns
% for a=linspace(0,0.2,ranger) %range of velocity
%     f=1; %adder rows
%     for l=linspace(10^-16,0.5,ranger)  %range of inter-fish distances dimensional (L) from 0 to 3BL
%         lam=l;
%         opts.Algorithm = 'levenberg-marquardt';
%         Root1 = fsolve(@SOLVEeq,X2,opts) %TSIDE BY SIDE
%         zita1 = Root1(1); %crea il vettore di root 1
%         th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
%         th2=Root1(3);
%         x = [zita1 -zita1 lam th1 th2];%variables to put in the jacobian to check stability and in the checkequations
%         E=equationmap(x)
%         EIGENMAP=jacob(x); % check jacobian with the values 
%         if E==1 && zita1>0.0001  %condition for side-by-side
%             R2(f,g) = Root1(1);
%             if EIGENMAP==1  %negative (STABLE SPIRAL)
%                 eigenva(f,g)=-1;
%             elseif EIGENMAP==2 %zero (MARGINALLY STABLE)
%                 eigenva(f,g)=0;
%             else %positive  (UNSTABLE)
%                 eigenva(f,g)=1;
%             end
%         else
%             eigenva(f,g)=NaN;
%             R2(f,g) = NaN;
%         end
%         f=f+1;
%     end
%     g=g+1;
% end
% a= linspace(0,0.2,ranger); % range of velocity (x axis)
% l= linspace((10^-16)/0.5,0.5/0.5,ranger);% range of inter-fish distance (from zero to 2BL) y axis
% [velocity,interfish]= meshgrid(a,l);
% cMap = [0, 0, 1; ...   % Blue for -1 positive 
%   0, 1, 1; ...       % Greenfor 0 zero
%   1, 0, 0]     % Red for 1 negative
% 
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
% %Show the height of the channel
% h=3.1/r  %mm because 3.1 is the dipole length in mm
% saveas(F6, 'staggeredK1R0.png')
% close(F6)