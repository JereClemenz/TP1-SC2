clc;clear all;close all;

%Defino las contantes del sistema
R=4.7e3;
L=10e-6;
C=100e-9;

%Defino las matrices
Mat_A=[-R/L -1/L;1/C 0];
Mat_B=[1/L;0];
Mat_C=[R 0;0 1];% Necesito el voltaje a la salida que es el producto entre la corriente de la malla y la resistencia

%Defino la entrada U
h=0.00001; %0.01 ms de paso
t = 0:h:0.01999; % Vector de tiempo de 0 a 20ms con paso de h=0.0001ms --> vector tamaño 2000
u = zeros(size(t)); % Vector de ceros de la misma longitud que el vector de tiempo
variable=0;%cuenta 1ms/h

for i = 1:length(t)
    
    if t(i) <= 0.005
        u(i)= 0; %primeros 5 ms vale cero
    elseif t(i) > 0.005 && variable<=(0.001/h)
        u(i) = 12;
        variable=variable+1;
    elseif t(i) > 0.005 && variable>(0.001/h)
        u(i) = -12;
        variable=variable+1;
    end
    
    if variable>(2*(0.001/h))
        variable=0;
    end
    
end
plot(t, u);
grid on;
title('u_t,V_a')
xlabel('Tiempo [Seg.]');
ylabel('Voltaje [Volt]');

%voltaje en el capacitor
Mat_C=[0 1];
sys1=ss(Mat_A,Mat_B,Mat_C,[]);
figure
lsim(sys1,u,t);
xlabel('Tiempo [Seg.]');
ylabel('Voltaje [Volt]');
title('Voltaje en el capacitor Vc(t)');

%corriente
Mat_C=[1 0];
sys2=ss(Mat_A,Mat_B,Mat_C,[]);
Corriente=lsim(sys2,u,t);
figure
plot(t,Corriente);
title('Corriente en el circuito i(t)');



