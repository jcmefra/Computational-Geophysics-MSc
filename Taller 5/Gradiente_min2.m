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
%Algoritmo gradiente-mínimos cuadrados
k=0;
m=zeros(256,1);
p=zeros(256,1);
beta=0;
s=-dw;
r=Gw'*s;
for k=1:1000
    p=-r+beta*p;
    alfa= (norm(r)*norm(r))/((p'*Gw')*(Gw*p));
    m=m+alfa*p;
    s=s+alfa*Gw*p;
    r_old=r;
    r=Gw'*s;
    beta=(norm(r)*norm(r))/(norm(r_old)*norm(r_old));
    res(k)=norm(((Gw*m)-dw));
    res1(k)=norm(m,1);
end
toc
figure, imagesc(reshape(m,16,16))
title('Matriz de slowness usando Gradiente Conjugado')
colormap bone
colorbar
figure, contourf(reshape(m,16,16))
title('Mapa de contornos de slowness usando Gradiente Conjugado')
colormap bone
colorbar
figure, semilogx(res)
title('Norma de los residuales por iteración')
xlabel('Iteraciones')
ylabel('Norma Residual ||Gm-d||_2')
figure, loglog(res1)
title('Norma de la solución por iteración')
xlabel('Iteraciones')
ylabel('Norma ||m||_2')

