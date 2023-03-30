clear all;clc;close all;

Da = readmatrix('Curvas_Medidas_Motor_2023.xlsx');%Da=datos

figure(1)
subplot(2,2,1)
plot(Da(:,1),Da(:,2));%grafico de velocidad angular
title('\omega_r')
xlabel('Tiempo [Seg.]');
ylabel('\omega [rad/s]');

subplot(2,2,2);
plot(Da(:,1),Da(:,3)); %grafico de corriente de armadura
title('i_a ');
xlabel('Tiempo [Seg.]');
ylabel('Corriente en armadura [A]');

subplot(2,2,3);
plot(Da(:,1),Da(:,4)); %grafico de tension de entrada
title('V_a ')
xlabel('Tiempo [Seg.]');
ylabel('V_a [V]');

subplot(2,2,4);
plot(Da(:,1),Da(:,5)); %grafico de corriente de armadura
title('T_L ')
xlabel('Tiempo [Seg.]');
ylabel('Torque [N*m]');

%Primero encontraremos Wr/Va
StepAmplitude=1; %ganancia encontrada de los datos de la curva RLC
K=198.2;
n=388;
dist=5000;
y_t1= Da(n,2);
t_t1=Da(n,1);
ii=1;
t_2t1=Da(n+dist,1);
y_2t1=Da(n+dist,2);
t_3t1=Da(n+2*dist,1);
y_3t1=Da(n+2*dist,2);

k1=(1/StepAmplitude)*y_t1/K-1;%Afecto el valor del Escalon
k2=(1/StepAmplitude)*y_2t1/K-1;
k3=(1/StepAmplitude)*y_3t1/K-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(be))/(sqrt(be));
%ingresar los t
T1_ang=-0.0005/log(alfa1);
T2_ang=-0.0005/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;

s=tf('s');
% G=(12*(T3_ang*s+1))/((T1_ang*s +1)*(T2_ang*s +1)) %le saco el cero, pq funciona asi
G1=(198.2)/((T1_ang*s +1)*(T2_ang*s +1))
%entrada de tension 12v
G1n=G1/12;%ft normalizada

% figure(5)
% hold on
% step(G1*exp(-s*0.025));
% plot(Da(:,1),Da(:,2),'r');
% title('encontrada');
% hold off

%Primero encontraremos Wr/TL -27 la velocidad angular
StepAmplitude=1; %ganancia encontrada de los datos de la curva RLC
Kgain2=171.4;% (400-(200+30)) aprox
n=15550;
dist=1100;
y_t1= Da(n,3);
t_t1=Da(n,1);
ii=1;
t_2t1=Da(n+dist,1);
y_2t1=Da(n+dist,3);
t_3t1=Da(n+2*dist,1);
y_3t1=Da(n+2*dist,3);

k1=(1/StepAmplitude)*(Kgain2-y_t1)/Kgain2-1;%Afecto el valor del Escalon
k2=(1/StepAmplitude)*(Kgain2-y_2t1)/Kgain2-1;
k3=(1/StepAmplitude)*(Kgain2-y_3t1)/Kgain2-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(be))/(sqrt(be));
%ingresar los t
T1_ang=-0.00005/log(alfa1);
T2_ang=-0.00005/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;

s=tf('s');
G2=(Kgain2*(T3_ang*s+1))/((T1_ang*s +1)*(T2_ang*s +1))
% G2=(Kgain2)/((T1_ang*s +1)*(T2_ang*s +1))
G2n=G2/(-1.04e-3)%Normalizo la FT TL de entrada -0.104Nm

figure(6)
hold on
step(143-G2*exp(-s*0.1501));
step(-198.2+G2*exp(-s*0.1501));
plot(Da(:,1),Da(:,2),'r');
title('encontrada wr/tl');
hold off

h=1e-5;
tsim=0.604;
t=0:h:(tsim-h);

%Defino de vuelta una entrada pq el paso es de 1e-7 y es mucha ccarga
%computacional

va=zeros(1,round(tsim/h));
for i=round(0.025/h):1:round(tsim/h)
    if i<=round(0.1501/h)
        va(1,i)=12;
    elseif i>round(0.1501/h)
        va(1,i)=-12;
    end
end
[y1 t1]=lsim(G1n,va(1,:),t);


%Defino el torque
tl=zeros(1,round(tsim/h));
for i=round(0.1504/h):1:round(tsim/h)
    tl(1,i)=-1.04e-3;
%     tl(1,i)=0;
end
[y2 t2]=lsim(G2n,tl(1,:),t);


figure(7)
hold on
plot(t,y1+y2);
plot(Da(:,1),Da(:,2),'--');
xlabel('Tiempo');
ylabel('wr');
title('Modelo dinámico encontrado');
legend({'Encontrada','Observada'},'Location','northwest','Orientation','horizontal')
hold off





%Corriente en funcion de la tension
% G3=G1*(0.00000033*s+0.0000001)
% figure(8)
% hold on
% step(G3*exp(-s*0.025),'b');
% plot(Da(:,1),Da(:,3),'r');
% title('encontrada ia/va');
% hold off
% %Normalizo
% G3n=G3/(12)
% 
% %Corriente en funcion del torque 
% zita=0.7834;
% tp=0.150122-0.1501002;
% wn=pi/(tp*sqrt(1-zita^2));
% G4=(0.105)*(wn^2)/(s^2+2*0.55*wn*s+wn^2)
% % figure(9)
% % hold on
% % step(G4*exp(-s*0.15),'b');
% % plot(Da(:,1),Da(:,3),'r');
% % title('encontrada ia/tl');
% % hold off
% %Normalizo
% G4n=G4/(0.104e-3)
% 
% [y3 t3]=lsim(G3n,va(1,:),t);
% [y4 t4]=lsim(G4n,tl(1,:),t);
% 
% figure(10)
% hold on
% plot(t,y3+y4);
% plot(Da(:,1),Da(:,3))
% xlabel('Tiempo');
% ylabel('ia');
% title('modelo dinámico encontrado para corriente. Real en naranja vs aproximado en azul');
% hold off
