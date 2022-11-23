close all
clear all
clc
%% test the new function
disp('=== FWI with LNSRCH === ');

% Inicialización de las características del modelo
Nx=210; 
Nz=68; 
offset_max = 4000; 
tEnd=2.5; 
dt=0.004;  t = 0:dt:tEnd-dt; frec=3;
Nt=tEnd/dt;  dz = 25; dx = dz;
vp_ori = zeros(Nx,Nz);
vp_ite = vp_ori;
Sx1=25;
Sx2=185;
Sx3=105;
Sx4=65;
Sx5=145;
Sz=3;
Sz1=3;
Sz2=3;
Sz3=3;
Sz4=3;
Sz5=3;
beta1=1; beta2=1; beta3=1; beta4=1; beta5=1;
f1=0; f2=0; f3=0; f4=0; f5=0;
grad1=zeros(Nx*Nz,1);
grad2=zeros(Nx*Nz,1);
grad3=zeros(Nx*Nz,1);
grad4=zeros(Nx*Nz,1);
grad5=zeros(Nx*Nz,1);

% Creación de la fuente
a = (pi*frec)^2; t0 = 1;
g1 = [-2*a*(t-t0).*exp(-a*(t-t0).^2)]';


%% 1) Creación del modelo de velocidades original
for iz=1:Nz; vp_ori(:,iz) = 2000+0.7*dz*(iz-1); end

figure
mesh(vp_ori)
view(-1,-1)

hold on

% Adquisición sobre el modelo original
%tic; [Pt_obs1]=propagator(vp_ori, g1, Sx1, Sz, dx, dz, dt, 4000, frec); toc % offset_max = 4000
%tic; [Pt_obs2]=propagator(vp_ori, g1, Sx2, Sz, dx, dz, dt, 4000, frec); toc % offset_max = 4000
tic; [Pt_obs3]=propagator(vp_ori, g1, Sx3, Sz, dx, dz, dt, 2125, frec); toc % offset_max = 2125
%tic; [Pt_obs4]=propagator(vp_ori, g1, Sx4, Sz, dx, dz, dt, 3125, frec); toc % offset_max = 3125
%tic; [Pt_obs5]=propagator(vp_ori, g1, Sx5, Sz, dx, dz, dt, 3125, frec); toc % offset_max = 3125

%% 2) Creación del modelo de velocidad inicial, para la primer iteración se
for iz=1:Nz; vp_ite(:,Nz-iz+1) = 3172.5-0.6*dz*(iz-1); end

mesh(vp_ite)
view(-1,-1)

for i=1:40
N(i) = getframe;    
vp_ite1=reshape(vp_ite,Nx*Nz,1);

%% 3) Cálculo del gradiente y función objetivo para la adquisición
% #1
%[f1,grad1] = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx1, Sz1, dx, dz, dt, 4000, frec, Pt_obs1, 1);% offset_max = 10000 k=6
% #2
%[f2,grad2] = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx2, Sz2, dx, dz, dt, 4000, frec, Pt_obs2, 2);% offset_max = 10000 k=6
% #3
[f3,grad3] = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx3, Sz3, dx, dz, dt, 2125, frec, Pt_obs3, 3);% offset_max = 10000 k=6
% #4
%[f4,grad4] = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx4, Sz4, dx, dz, dt, 3125, frec, Pt_obs4, 4);% offset_max = 10000 k=6
% #5
%[f5,grad5] = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx5, Sz5, dx, dz, dt, 3125, frec, Pt_obs5, 5);% offset_max = 10000 k=6


%% 4) Función objetivo total y gradiente total
f(i) = f1 + f2 + f3 + f4 + f5

%% 5) Aplicando el factor de escala
grad=beta1*grad1 + beta2*grad2 + beta3*grad3 + beta4*grad4 + beta5*grad5;    

%% 6) Actualizando el modelo de velocidades
vp_ite = vp_ite + reshape(grad(:,:),Nx,Nz);
        
min(min(vp_ite))
mesh(vp_ite)
view(-1,-1)
end
movie(N,4)