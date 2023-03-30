clear all;clc;close all;

Da = readmatrix('Curvas_Medidas_RLC.xls');%Da=datos

figure(1)
plot(Da(:,1),Da(:,2)); %grafico de corriente
title('i_a variable de estado x_1')
xlabel('Tiempo [Seg.]');
ylabel('Corriente [Amp.]');

figure(2)
plot(Da(:,1),Da(:,3)); %grafico de voltaje de capacitor
title('V_c var. de est. x_2')
grid on;
xlabel('Tiempo [Seg.]');
ylabel('Voltaje [Volt]');


StepAmplitude=1; %ganancia encontrada de los datos de la curva RLC
K=12;
n=151;
dist=50;
y_t1= Da(n,3);
t_t1=Da(n,1);
ii=1;

t_2t1=Da(n+dist,1);
y_2t1=Da(n+dist,3);

t_3t1=Da(n+2*dist,1);
y_3t1=Da(n+2*dist,3);

k1=(1/StepAmplitude)*y_t1/K-1;%Afecto el valor del Escalon
k2=(1/StepAmplitude)*y_2t1/K-1;
k3=(1/StepAmplitude)*y_3t1/K-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(k1+alfa2)/(alfa1-alfa2);

T1_ang=-0.005/log(alfa1);
T2_ang=-0.005/log(alfa2);
% T1_ang=-y_t1/log(alfa1);
% T2_ang=-y_t1/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;


s=tf('s');
%G=(12*(T3_ang*s+1))/((T1_ang*s +1)*(T2_ang*s +1)) %le saco el cero, pq funciona asi
G=(12)/((T1_ang*s +1)*(T2_ang*s +1))
[numG,denG]=tfdata(G,'v');

figure(3)
hold on
step(G*exp(-s*0.01),'y');
plot(Da(:,1),Da(:,3),'--r');
title('encontrada');
legend({'Encontrada','Observada'},'Location','northwest','Orientation','horizontal')
hold off

%Encontramos los parametros RLC
%Para este caso iteraremos el valor de R ya que el primer coeficiente de la
%ft encontrada es LC y el segundo 2RC.
R=220;
Cap=denG(1,2)/(R);
L=denG(1,1)/Cap;


% Armamos la tension de entrada u, que para este caso observando los datos
% aportados tenemos una entrada: 10ms 0v 40ms 12v 50ms -12v
ts=0.1;%tiempo de simulacion
h=0.00001; % paso entre muestras 0.01ms
t=0:h:(ts-h);% vector tiempo

u=zeros(1,ts/h);

for i=round(0.01/h):1:round(ts/h)
    
    if i>(0.01/h+1) && i<=(0.05/h)
        u(1,i)=12;
    elseif i>(0.05/h) && i<=(0.1/h)
        u(1,i)=-12;
    end
end
% figure(4)
% plot(t,u); %grafico de entrada u
% title('Entrada V')
% xlabel('Tiempo [Seg.]');
% ylabel('Voltaje [Volt]');

%Definimos las matrices
A=[-R/L -1/L; 1/Cap 0];
B=[1/L; 0];
C=[1 0];

%Definimos el sistema para corriente
sys1=ss(A,B,C,[]);
y_i=lsim(sys1,u,t);
figure(5)
hold on
grid on;
xlabel('Tiempo [Seg.]');
ylabel('Corriente [Amp.]');
plot(t,y_i,'.-r');
plot(Da(:,1),Da(:,2),'g');
legend({'Encontrada','Observada'},'Location','northwest','Orientation','horizontal')
hold off

%Definimos el sistema para tension
C=[0 1]
sys2=ss(A,B,C,[]);
y_vc=lsim(sys2,u,t);
figure(6)
hold on
grid on;
xlabel('Tiempo [Seg.]');
ylabel('Voltaje [Volt]');
plot(t,y_vc,'.-r');
plot(Da(:,1),Da(:,3),'g');
legend({'Encontrada','Observada'},'Location','northwest','Orientation','horizontal')
hold off
