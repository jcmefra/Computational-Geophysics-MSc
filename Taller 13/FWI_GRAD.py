
import numpy as np
from numpy import linalg as LA

from propagator import propagator
from propagator1 import propagator1

def FWI_GRAD(x, Nx, Nz, Nt, g1, Sx1, Sz, dx, dz, dt, offset_max, frec, Pt_obs,k):

    # Propagación del frente de onda sobre el modelo de velocidad
    Pt_mod,P_mod,d2P_dt2 = propagator(np.reshape(x[:],(Nx,Nz),order='F'), g1, Sx1, Sz, dx, dz, dt, offset_max, frec)


    # 1) Cálculo de la función en el punto
    # tomando como referencia p5.pdf
    if k < 3:
        ad = 5
    elif k == 3 or k == 5:
        ad = 1
    else:
        ad = 0 

    
    f = np.dot(0.5*(np.reshape(Pt_mod,(Nt*(Nx-40-ad),1),order='F') - np.reshape(Pt_obs,(Nt*(Nx-40-ad),1),order='F')).T,(np.reshape(Pt_mod,(Nt*(Nx-40-ad),1),order='F') - np.reshape(Pt_obs,(Nt*(Nx-40-ad),1),order='F')))
    
    # 2) Cálculo y visualización del residual
    res = Pt_mod-Pt_obs


    # 3) Retropropagación del residual
    Pt_back = propagator1(np.reshape(x[:],(Nx,Nz),order='F'), np.flipud(res), np.arange(20,Nx-60,1), Sz, dx, dz, dt,offset_max, frec)


    P_back_t = np.zeros((Nx,Nz,Nt))

    for it in range(0,Nt-1):
        P_back_t[:,:,it]=Pt_back[1][:,:,Nt-it-1]
    

    # 5) Cálculo del gradiente normal
    gradient = -dt*np.sum(P_back_t*d2P_dt2,axis=2)

    gradient[Sx1-1,Sz-1]=gradient[Sx1-2,Sz-1]
    gradient[Sx1-1,Sz]=gradient[Sx1-2,Sz-1]
    gradient=gradient*(1/LA.norm(gradient,2))*100

    g=np.reshape(gradient,(Nx*Nz,1),order='F')

    # Recomendación profesora Clara
    
    # Visualización del gradiente
    #  figure; 
    #  imagesc(1e-3*dx*(1:Nx),1e-3*dz*(1:Nz),gradient')
    #  caxis(5*[-1 1]), axis equal, axis tight, grid on, axis([0 5.2 0 1.7]), colorbar
    #  xlabel('Distance (km)'), ylabel('Depth (km)'), title('Gradient'), drawnow
    
    #  figure
    #  mesh(gradient)
    #  view(0,0)
    return f,g