function [f,g] = FWI_GRAD(x, Nx, Nz, Nt, g1, Sx1, Sz, dx, dz, dt, offset_max, frec, Pt_obs,k)

% Propagación del frente de onda sobre el modelo de velocidad
tic; [Pt_mod,P_mod,d2P_dt2]=propagator(reshape(x(:,:),Nx,Nz), g1, Sx1, Sz, dx, dz, dt, offset_max, frec); toc


% 1) Cálculo de la función en el punto
% tomando como referencia p5.pdf
    if k<3
        ad=5;
    elseif k==3 || k==5
        ad=1;
    else
        ad=0;
    end
    
    f= 0.5*(reshape(Pt_mod,Nt*(Nx-40-ad),1) - reshape(Pt_obs,Nt*(Nx-40-ad),1))'*(reshape(Pt_mod,Nt*(Nx-40-ad),1) - reshape(Pt_obs,Nt*(Nx-40-ad),1));

% 2) Cálculo y visualización del residual
res =Pt_mod-Pt_obs;


% 3) Retropropagación del residual
tic; [Pt_back,P_back]=propagator(reshape(x(:,:),Nx,Nz), flipud(res), 21:Nx-60, Sz, dx, dz, dt,offset_max, frec); toc


P_back_t = zeros(Nx,Nz,Nt);

for it=1:Nt; P_back_t(:,:,it)=P_back(:,:,Nt-it+1); end

%% 5) Cálculo del gradiente normal
gradient = -dt*sum(P_back_t.*d2P_dt2,3);

gradient(Sx1,Sz)=gradient(Sx1-1,Sz);
gradient(Sx1,Sz+1)=gradient(Sx1-1,Sz);
gradient=gradient*(1/norm(gradient))*100;

g=reshape(gradient,Nx*Nz,1);

%% Recomendación profesora Clara

% Visualización del gradiente
%  figure; 
%  imagesc(1e-3*dx*(1:Nx),1e-3*dz*(1:Nz),gradient')
%  caxis(5*[-1 1]), axis equal, axis tight, grid on, axis([0 5.2 0 1.7]), colorbar
%  xlabel('Distance (km)'), ylabel('Depth (km)'), title('Gradient'), drawnow

%  figure
%  mesh(gradient)
%  view(0,0)
