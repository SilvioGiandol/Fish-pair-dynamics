%CROSS STREAM EQUILIBRIA SOLVED WITH F SOLVE
% Author:  S. Giandola, 05/31/2022
clear all, clc, close all
% 
%% START PLOTTING simple
global   r kappa a l kappap kappav ps kappas

kappas=0;
kappap=0; % VISION
kappav=0; % VISION
kappa=0;  %NO LATERAL-LINE SENSING

l=0.9; % inter-fish distance
r=0.1;
% Initial guess for numerical approach
ps=1;

in1=0.33;
in2=-0.33;
in3=0;
in4=0.1;
in5=-0.1; 
X1=[in1 pi pi]; 
X2=[in2 pi pi]; 
X3=[in3 pi pi]; 
X4=[in4 pi pi]; 
X5=[in5 pi pi]; 

disp('CROSS STREAM EQUILIBRIA')

af=0.16;  %velocità finale, per check completo metterei 0.5

a= 0:0.01:af;  
R1 = zeros(length(a),1);  %root 1
R2 = zeros(length(a),1); %root 2
R3 = zeros(length(a),1); %root 3
R4 = zeros(length(a),1); %root 4
R5 = zeros(length(a),1); %root 5

Red1 = zeros(length(a),1); %unstable  
Red2 = zeros(length(a),1); %unstable
Red3 = zeros(length(a),1); %unstable

velocity= zeros(length(a),1); % Non dimensional velocity: alpha vector

Th1root1=zeros(length(a),1); % Per vedere theta1
Th2root1=zeros(length(a),1); % Per vedere theta1

Th1root2=zeros(length(a),1); % Per vedere theta1
Th2root2=zeros(length(a),1); % Per vedere theta1

Th1root4=zeros(length(a),1); % Per vedere theta1
Th2root4=zeros(length(a),1); % Per vedere theta1

Th1root5=zeros(length(a),1); % Per vedere theta1
Th2root5=zeros(length(a),1); % Per vedere theta1

Th1root3=zeros(length(a),1); % Per vedere theta1
Th2root3=zeros(length(a),1); % Per vedere theta1
f=1;

for a=0:0.01:af  %range of alpha
    lam=l;
    %Create arrays of solution to plot:
    opts.Algorithm = 'levenberg-marquardt';
    %ROOT1:
    Root1 = fsolve(@SOLVEeq,X1,opts);
    zita1 = Root1(1); %crea il vettore di root 1
    th1=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
    th2=Root1(3);
    
    x=[zita1 -zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
    J=jacob(x) %check jacobian with the values 
    E=equationmap(x) %check if equations are satisfied
    %if J==2 && E==1 
    if E==1  %of it's a fixed points and it is also stable
        R1(f) = Root1(1); %crea il vettore di root 1
        Th1root1(f)=Root1(2);%siccome theta potrebbe non essere pi è meglio checkare
        Th2root1(f)=Root1(3);
    else
        R1(f) = NaN; %crea il vettore di root 1
        Th1root1(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
        Th2root1(f)=NaN;
    end
    %ROOT2:
    Root2 = fsolve(@SOLVEeq,X2,opts);
    zita1 = Root2(1); %crea il vettore di root 1
    th1=Root2(2);%siccome theta potrebbe non essere pi è meglio checkare
    th2=Root2(3);
    x=[zita1 -zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
    J=jacob(x) %check jacobian with the values 
    E=equationmap(x) %check if equations are satisfied
    %if J==2 && E==1 
    if E==1  %of it's a fixed points and it is also stable
        R2(f) = Root2(1); %crea il vettore di root 2
        Th1root2(f)=Root2(2);%siccome theta potrebbe non essere pi è meglio checkare
        Th2root2(f)=Root2(3);
    else
        R2(f) = NaN; %crea il vettore di root 1
        Th1root2(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
        Th2root2(f)=NaN;
    end

    %ROOT4:
    Root4 = fsolve(@SOLVEeq,X4,opts)
    zita1 = Root4(1); %crea il vettore di root 1
    th1=Root4(2);%siccome theta potrebbe non essere pi è meglio checkare
    th2=Root4(3);
    x=[zita1 -zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
    J=jacob(x) %check jacobian with the values 
    E=equationmap(x) %check if equations are satisfied
    if J==2 && E==1 
    %if E==1  %of it's a fixed points and it is also stable
        R4(f) = Root4(1); %crea il vettore di root 5
        Th1root4(f)=Root4(2);%siccome theta potrebbe non essere pi è meglio checkare
        Th2root4(f)=Root4(3);
        Red1(f)=NaN;
    else
        R4(f) = NaN; %crea il vettore di root 1
        Red1(f)=Root4(1);
        Th1root4(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
        Th2root4(f)=NaN;
    end
    %ROOT5:
    Root5 = fsolve(@SOLVEeq,X5,opts)
    zita1 = Root5(1); %crea il vettore di root 1
    th1=Root5(2);%siccome theta potrebbe non essere pi è meglio checkare
    th2=Root5(3);
    x=[zita1 -zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
    J=jacob(x) %check jacobian with the values 
    E=equationmap(x) %check if equations are satisfied
    if J==2 && E==1 
    %if E==1  %of it's a fixed points and it is also stable
        R5(f) = Root5(1); %crea il vettore di root 5
        Th1root5(f)=Root5(2);%siccome theta potrebbe non essere pi è meglio checkare
        Th2root5(f)=Root5(3);
        Red2(f) = NaN;
       
    else
        R5(f) = NaN; %crea il vettore di root 1
        Red2(f)=Root5(1);
        Th1root5(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
        Th2root5(f)=NaN;
    end

    %ROOT3 (CLOSE TO ZERO)
    Root3 = fsolve(@SOLVEeq,X3,opts);
    zita1 = Root3(1); %crea il vettore di root 1
    th1=Root3(2);%siccome theta potrebbe non essere pi è meglio checkare
    th2=Root3(3);
    x=[zita1 -zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
    J=jacob(x) %check jacobian with the values 
    E=equationmap(x) %check if equations are satisfied
    if J==2 && E==1 
       R3(f) = Root3(1);
       Th1root3(f)=Root3(2);
       Th2root3(f)=Root3(3);
       Red3(f) = NaN;
    else
       Red3(f) = Root3(1);
       Th1root3(f)=NaN;
       Th2root3(f)=NaN;
    end
    velocity(f)= a;
    a
    f=f+1;
end

figure('DefaultAxesFontSize',16)
set(1,'units','normalized','name','Cross-stream equilibria','numbertitle','off');

plot(velocity,R3,'g','LineWidth',3) %da rendere rosso quando arriviamo a lambda 0.4

hold on
plot(velocity,Red1,'r',velocity,Red2,'r',velocity,Red3,'r','LineWidth',2.5,'LineWidth',3)
plot(velocity,R1,'r',velocity,R2,'r',velocity,R3,'r','LineWidth',2.5,'LineWidth',3)
%plot(velocity,Red2,'r','LineWidth',2.5,'LineWidth',3)
plot(velocity,R5,'g',velocity,R4,'g','LineWidth',3)
%plot(velocity,R1,'r',velocity,R2,'r','LineWidth',3)

axis([0 af -0.5 0.5])

xlabel('\alpha','FontSize', 20),ylabel('\xi','FontSize', 20)

title('Cross stream equilibria ')

%TABLE TO SHOW THE VARIATION
%data = table(velocity,R1,R2,R3,ratio3,R4,Th1root4,Th2root4,ratio4,R5,ratio5,Th1root5,Th2root5)
data = table(velocity,R3,Red2) % to show the alpha critical
%data2=table(velocity,R4,Th1root4,Th2root4,ratio4)

% %% START PLOTTING STAGGERED WITHOUT CONSTRAINTS on ksi GENERALIZED
% 
% global r kappa a l kappap kappav ps kappas
% kappas=0;
% kappap=0; % W/OUT VISION
% kappav=0; % W/OUT VISION
% kappa=0; % NO LATERAL-LINE SENSING
% l=0.8;
% r=0.1; 
% %Initial guess for numerical approach
% p=1;
% 
% in1=0.33;
% in2=-0.33;
% in3=0;
% in4=0.09;
% in5=-0.09;  %0.07 for lambda=0.1 e 0.2 mentre 0.18 for the other lambda
% 
% X1=[in1 in2 pi pi]; 
% X2=[in2 in1 pi pi]; 
% X3=[in3 in3 pi pi]; 
% X4=[in4 in5 pi pi]; 
% X5=[in5 in4 pi pi]; 
% 
% disp('CROSS STREAM EQUILIBRIA')
% 
% af=0.16; %velocità finale, per check completo metterei 0.5
% 
% % a= 0:0.01:af; 
% a= 0:0.01:af;
% ZITA1R1 = zeros(length(a),1);  %root 1
% ZITA2R1 = zeros(length(a),1);  %root 1
% ZITA1R2 = zeros(length(a),1);  %root 1
% ZITA2R2 = zeros(length(a),1);  %root 1
% ZITA1R3 = zeros(length(a),1);  %root 1
% ZITA2R3= zeros(length(a),1);  %root 1
% ZITA1R4 = zeros(length(a),1);  %root 1
% ZITA2R4= zeros(length(a),1);  %root 1
% ZITA1R5= zeros(length(a),1);  %root 1
% ZITA2R5 = zeros(length(a),1);  %root 1
% 
% velocity= zeros(length(a),1); % Non dimensional velocity: alpha vector
% 
% Th1root1=zeros(length(a),1); % Per vedere theta1
% Th2root1=zeros(length(a),1); % Per vedere theta1
% ratio1=zeros(length(a),1);
% 
% Th1root2=zeros(length(a),1); % Per vedere theta1
% Th2root2=zeros(length(a),1); % Per vedere theta1
% ratio2=zeros(length(a),1);
% 
% Th1root4=zeros(length(a),1); % Per vedere theta1
% Th2root4=zeros(length(a),1); % Per vedere theta1
% ratio4=zeros(length(a),1);
% 
% Th1root5=zeros(length(a),1); % Per vedere theta1
% Th2root5=zeros(length(a),1); % Per vedere theta1
% ratio5=zeros(length(a),1);
% 
% Th1root3=zeros(length(a),1); % Per vedere theta1
% Th2root3=zeros(length(a),1); % Per vedere theta1
% ratio3=zeros(length(a),1);
% 
% f=1;
% 
% for a=0:0.01:af %range of alpha
%     lam=l;
%     ps=p;
%     %Create arrays of solution to plot:
%     opts.Algorithm = 'levenberg-marquardt';
%     %ROOT1:
%     Root1 = fsolve(@SOLVEeq,X1,opts)
%     zita1 = Root1(1); %crea il vettore di root 1
%     zita2=Root1(2);
%     th1=Root1(3);%siccome theta potrebbe non essere pi è meglio checkare
%     th2=Root1(4);
%     
%     x=[zita1 zita2 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
%     J=jacob(x) %check jacobian with the values 
%     E=equationmap(x) %check if equations are satisfied
%     %E=1;
%     %if J==1 && E==1 
%     if E==1  %of it's a fixed points and it is also stable
%         ZITA1R1(f) = Root1(1); %crea il vettore di root 1
%         ZITA2R1(f) = Root1(2);
%         Th1root1(f)=Root1(3);%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root1(f)=Root1(4);
%         ratio1(f)=p;
%     else
%         ZITA1R1(f) = NaN; %crea il vettore di root 1
%         ZITA2R1(f) = NaN;
%         Th1root1(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root1(f)=NaN;
%         ratio1(f)=NaN;
%     end
%     %ROOT2:
%     Root2 = fsolve(@SOLVEeq,X2,opts)
%     zita1 = Root2(1); %crea il vettore di root 1
%     zita2= Root2(2);
%     th1=Root2(3);%siccome theta potrebbe non essere pi è meglio checkare
%     th2=Root2(4);
%     x=[zita1 zita2 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
%     J=jacob(x) %check jacobian with the values 
%     E=equationmap(x) %check if equations are satisfied
%     %E=1;
%     %if J==1 && E==1 
%     if E==1  %of it's a fixed points and it is also stable
%         ZITA1R2(f) = Root2(1); %crea il vettore di root 1
%         ZITA2R2(f) = Root2(2);
%         Th1root2(f)=Root2(3);%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root2(f)=Root2(4);
%         ratio2(f)=p;
%     else
%         ZITA1R2(f) = NaN; %crea il vettore di root 1
%         ZITA2R2(f) = NaN;
%         Th1root2(f) = NaN;%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root2(f) = NaN;
%         ratio2(f) = NaN;
%     end
% 
%     %ROOT4:
%     Root4 = fsolve(@SOLVEeq,X4,opts)
%     zita1 = Root4(1); %crea il vettore di root 1
%     zita2=Root4(2);
%     th1=Root4(3);%siccome theta potrebbe non essere pi è meglio checkare
%     th2=Root4(4);
%     x=[zita1 zita2 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
%     J=jacob(x) %check jacobian with the values 
%     E=equationmap(x) %check if equations are satisfied
%     %E=1;
%     if J==1 && E==1 
%     %if E==1  %of it's a fixed points and it is also stable
%         ZITA1R4(f) = Root4(1); %crea il vettore di root 1
%         ZITA2R4(f) = Root4(2);
%         Th1root4(f)=Root4(3);%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root4(f)=Root4(4);
%         ratio4(f)=p;
%     else
%         ZITA1R4(f) = NaN; %crea il vettore di root 1
%         ZITA2R4(f) = NaN;
%         Th1root4(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root4(f)=NaN;
%         ratio4(f)=NaN;
%     end
% 
%     %ROOT5:
%     Root5 = fsolve(@SOLVEeq,X5,opts)
%     zita1 = Root5(1); %crea il vettore di root 1
%     zita2=Root5(2);
%     th1=Root5(3);%siccome theta potrebbe non essere pi è meglio checkare
%     th2=Root5(4);
%     x=[zita1 zita2 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
%     J=jacob(x) %check jacobian with the values 
%     E=equationmap(x) %check if equations are satisfied
%     %E=1;
%     if J==1 && E==1 
%     %if E==1  %of it's a fixed points and it is also stable
%         ZITA1R5(f) = Root5(1); %crea il vettore di root 1
%         ZITA2R5(f) = Root5(2);
%         Th1root5(f)=Root5(3); %siccome theta potrebbe non essere pi è meglio checkare
%         Th2root5(f)=Root5(4);
%         ratio5(f)=p;
%     else
%         ZITA1R5(f) = NaN; %crea il vettore di root 1
%         ZITA2R5(f) = NaN;
%         Th1root5(f) = NaN;%siccome theta potrebbe non essere pi è meglio checkare
%         Th2root5(f) = NaN;
%         ratio5(f) = NaN;
%     end
% 
% %    %ROOT3 (CLOSE TO ZERO)
% %     Root3 = fsolve(@SOLVEeq,X3,opts)
% %     zita1 = real(Root3(1)); %crea il vettore di root 1
% %     zita2=real(Root3(2));
% %     th1=real(Root3(3));%siccome theta potrebbe non essere pi è meglio checkare
% %     th2=real(Root3(4));
% %     x=[zita1 zita1 lam th1 th2]; %variables to put in the jacobian to check stability and in the checkequations
% %     J=jacob(x) %check jacobian with the values 
% %     E=equationmap(x) %check if equations are satisfied
% %     %E=1;
% %     if J==1 && E==1 
% %     %if E==1  %of it's a fixed points and it is also stable
% %         ZITA1R3(f) = Root3(1); %crea il vettore di root 1
% %         ZITA2R3(f) = Root3(2);
% %         Th1root3(f)=Root3(3);%siccome theta potrebbe non essere pi è meglio checkare
% %         Th2root3(f)=Root3(4);
% %         ratio3(f)=p;
% %     else
% %         ZITA1R3(f) = NaN; %crea il vettore di root 1
% %         ZITA2R3(f) = NaN;
% %         Th1root3(f)=NaN;%siccome theta potrebbe non essere pi è meglio checkare
% %         Th2root3(f)=NaN;
% %         ratio3(f)=NaN;
% %     end
%     velocity(f)= a;
%     a
%     f=f+1;
% end
% 
% 
% figure(1)
% set(1,'units','normalized','name','Cross-stream equilibria 1','numbertitle','off');
% % 
% % plot(velocity,ZITA1R3,'g','LineWidth',3) %da rendere rosso quando arriviamo a lambda 0.4
% 
% hold on
% 
% % plot(velocity,R1,'r',velocity,R2,'r',velocity,R3,'r','LineWidth',2.5,'LineWidth',3)
% 
% plot(velocity,ZITA1R5,'g',velocity,ZITA1R4,'g','LineWidth',3)
% plot(velocity,ZITA1R1,'r',velocity,ZITA1R2,'r','LineWidth',3)
% 
% axis([0 af -0.5 0.5])
% 
% xlabel('\alpha','FontSize', 20),ylabel('\xi_1','FontSize', 20)
% 
% title('CROSS-STREAM EQUILIBRIA')
% figure(2)
% set(2,'units','normalized','name','Cross-stream equilibria 2','numbertitle','off')
% plot(velocity,ZITA2R3,'g','LineWidth',3) %da rendere rosso quando arriviamo a lambda 0.4
% 
% hold on
% 
% % plot(velocity,R1,'r',velocity,R2,'r',velocity,R3,'r','LineWidth',2.5,'LineWidth',3)
% 
% plot(velocity,ZITA2R5,'g',velocity,ZITA2R4,'g','LineWidth',3)
% plot(velocity,ZITA2R1,'r',velocity,ZITA2R2,'r','LineWidth',3)
% 
% axis([0 af -0.5 0.5])
% 
% xlabel('\alpha','FontSize', 20),ylabel('\xi_2','FontSize', 20)
% title('CROSS-STREAM EQUILIBRIA')
% 
% %TABLE TO SHOW THE VARIATION
% 
% data = table(velocity,ZITA1R1,ZITA1R2,ZITA1R3,ratio3,ZITA1R4,Th1root4,Th2root4,ratio4,ZITA1R5,ratio5,Th1root5,Th2root5)
% data2=table(velocity,ZITA1R4,Th1root4,Th2root4,ratio4)
% data3=table(velocity,ZITA1R3)

