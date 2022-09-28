% Taller 6
clc
clear all
close all
%Valores medidos
hh=[0.72 0.49 0.30 0.20 0.16 0.12];
tt=[5.0 10.0 20.0 30.0 40.0 50.0];
sigma=0.01;
%Parametros medidos
Q=50;
d=60;
%Ecuacion no lineal
syms S T t;
h=(Q/(4*pi*T*t))*exp((-(d^2)*S)/(4*T*t));
%Derivadas
dS=diff(h,S);
dT=diff(h,T);
%Jacobiano
n=length (hh);
tic
for i=1:n
    J(i,1)=subs(dS,t,tt(i))/sigma;
    J(i,2)=subs(dT,t,tt(i))/sigma;
end
%F_vector
for i=1:n
    F(i)=(subs(h,t,tt(i))-hh(i))/sigma;
end
%Condiciones iniciales[S,T].
S=0.005;   %0.003;
T=0.5; %0.5;
umbral=1e-06;
Deltam=[1 1];
iter = 0;
while Deltam(1)>umbral || Deltam(2)>umbral
    iter = iter+1;
    %Evalua J y F
    J1=eval(J);
    F1=eval(F);
    %3y4
    A=J1'*J1;
    b=-J1'*F1';
    Ab=[A b];
    sol=rref(Ab);
    Deltam=sol(:,3);
    S=S+Deltam(1);
    T=T+Deltam(2);
end
% Aspectos Estadisticos
cov=inv(A);
incer=-1.96*sqrt(diag(cov));
S1=S+incer(1);
S2=S-incer(1);
inter_S=[S1 S2];
T1=T+incer(2);
T2=T-incer(2);
inter_T=[T1 T2];
for i=1:n
    m(i)=eval(subs(h,t,tt(i)));
end
toc
%elipsoide de error
Vector_m = [S;T];
Delta2 = chi2inv(0.95,2);
plot_ellipse(Delta2,cov,Vector_m,10000)
% Grafica
figure, plot(tt, hh,'b*')
hold on
grid on
plot(tt,m,'ro')
legend('obs', 'mod')
xlabel('Tiempo(h)')
ylabel('Elevación de la tabla de agua(m)')
title("Estimados vs Observaciones")
