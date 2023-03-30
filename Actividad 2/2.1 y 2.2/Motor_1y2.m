clear all;close all;clc;

%Definicion de parametros

Laa=366e-6;
J=5e-9;
Ra=55.6;
B=0;
Ki=6.49e-3;
Km=6.53e-3;

A=[-Ra/Laa -Km/Laa 0;Ki/J -B/J 0; 0 1 0];
B=[1/Laa 0;0 -1/J;0 0];

%Definiion de vectores y Dt=h
tsim=5;
h=1e-5;
t=0:h:(tsim-h);

%Definicion de entrada U1=Va
u1=zeros(1,round(tsim/h));
for i=(round(0.5/h)):1:(tsim/h+1)
    u1(1,i)=12;
end
figure(1)
plot(t,u1);
title('Entrada U_1=V_a de 12V');

%Definicion de entrada de TL=U2 u2(1,i)=0.00140072;
u2=zeros(1,round(tsim/h));
for i=(3/h):1:(tsim/h+1)
    u2(1,i)=0.00140072;
%       u2(1,i)=0;
end
figure(2)
plot(t,u2);
title('Entrada U_2=T_L');

%Simulacion 
%condiciones iniciales
x(1,1)=0;
x(2,1)=0;
x(3,1)=0;
u=[u1;u2];

% B1=[1/Laa;0;0];
% B2=[0;-1/J;0];
for i=1:1:(tsim/h+1)
    %Variables del sistema lineal
    x1(1,i)= x(1,1);
    x2(1,i)= x(2,1);
    x3(1,i)= x(3,1);
    %Sistema lineal
    xp=A*x+B*u(:,i);
%   xp=A*x+B1*u1(1,i)+B2*u2(1,i);
    x=x+h*xp;
    
end


figure(3)
plot(t,x1);
title('Corriente de armadura i_a');
xlabel('Tiempo (seg.)');
ylabel('Corriente (A)');
grid on;

figure(4)
plot(t,x2);
title('Velocidad angular \omega_r');
xlabel('Tiempo (seg.)');
ylabel('Velocidad angular (rad/s)');
grid on;

figure(5)
plot(t,x3);
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;











