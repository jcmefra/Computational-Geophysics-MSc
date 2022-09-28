close all
clear all
clc
load ('jitorres_crosswell.mat')
sigma=0.5;

tic
for i=1:256
v(i)=1/sigma ;
end
%Matriz
W=diag(v,0);
Gw=W*G;
dw=W*dn;
%Algoritmo de Kaczmarz
m=zeros(256,1);
residual=1;
%contador=0;
%while residual>0.0001
for j=1:5221
for i=1:255
    m=m-((((Gw(i,:)*m)-dw(i))/(Gw(i,:)*Gw(i,:)'))*(Gw(i,:)'));
end
%contador=contador+1;
%residual=norm(((Gw*m)-dw));
residual(j)=norm(((Gw*m)-dw));
res1(j)=norm(m,1);
end
toc

figure, imagesc(reshape(m,16,16))
title('Matriz de slowness usando Kaczmarz')
colormap bone
colorbar
figure, contourf(reshape(m,16,16))
title('Mapa de contornos de slowness usando Kaczmarz')
colormap bone
colorbar
figure, semilogx(residual)
title('Norma de los residuales por iteración')
xlabel('Iteraciones')
ylabel('Norma Residual ||Gm-d||_2')
figure, loglog(res1)
title('Norma de la solución por iteración')
xlabel('Iteraciones')
ylabel('Norma ||m||_2')


