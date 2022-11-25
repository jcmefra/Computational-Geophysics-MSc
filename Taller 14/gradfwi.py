from FWI_GRAD2 import FWI_GRAD
from propagator import propagator
import numpy as np

Nx = 210
Nz = 68  
tEnd = 2.5 
dt = 0.004
Nt = tEnd/dt 
Nt= int(Nt)
Sx1 = 25
Sx2 = 185
Sx3 = 105
Sx4 = 65
Sx5 = 145
Sz1 = 3
Sz2 = 3
Sz3 = 3
Sz4 = 3
Sz5 = 3
Sz = 3
dz = 25
dx = 25
t = np.arange(0,tEnd,dt)
frec = 3

a = (np.pi*frec)**2
t0 = 1
g1 = (-2*a*(t-t0)*np.exp(-a*(t-t0)**2)).T
g1 = np.reshape(g1,(np.size(g1),1))

vp_ori = np.zeros((Nx,Nz))

for iz in range(0,Nz):
    vp_ori[:,iz] = 2000+0.7*dz*(iz)

def gradfwi(vp_ite1):
    Pt_obs1 = propagator(vp_ori, g1, Sx1, Sz, dx, dz, dt, 4000, frec)
    
    Pt_obs2 = propagator(vp_ori, g1, Sx2, Sz, dx, dz, dt, 4000, frec)
    
    Pt_obs3 = propagator(vp_ori, g1, Sx3, Sz, dx, dz, dt, 2125, frec)
    
    Pt_obs4 = propagator(vp_ori, g1, Sx4, Sz, dx, dz, dt, 3125, frec)
    
    Pt_obs5 = propagator(vp_ori, g1, Sx5, Sz, dx, dz, dt, 3125, frec)
    
    #vp_ite1 = np.reshape(vp_ite,(Nx*Nz), order='F')
    
    f1,grad1 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx1, Sz1, dx, dz, dt, 4000, frec, Pt_obs1[0], 1) # offset_max = 10000 k=6
    #2
    f2,grad2 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx2, Sz2, dx, dz, dt, 4000, frec, Pt_obs2[0], 2) # offset_max = 10000 k=6
    ##3
    f3,grad3 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx3, Sz3, dx, dz, dt, 2125, frec, Pt_obs3[0], 3) # offset_max = 10000 k=6
    #4
    f4,grad4 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx4, Sz4, dx, dz, dt, 3125, frec, Pt_obs4[0], 4) # offset_max = 10000 k=6
    #5
    f5,grad5 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx5, Sz5, dx, dz, dt, 3125, frec, Pt_obs5[0], 5) # offset_max = 10000 k=6

    grad = grad1 + grad2 + grad3 + grad4 + grad5 
    
    # vuelvo matriz para eliminar puntos de mayor energía ocasionados por las fuentes
    grad1 = np.reshape(grad[:],(Nx,Nz),order='F')
    # pixeles hacia abajo a partir de la posicion de la fuente
    offset= 0 # 0,1,2,3,4
    grad1[0:Sz1+1+offset,:]=0
    
    # vuelvo vector para retornar salida
    grad = np.reshape(grad1,(Nx*Nz), order='F')
    return grad