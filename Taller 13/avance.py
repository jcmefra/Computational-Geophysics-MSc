import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from propagator import propagator
from FWI_GRAD import FWI_GRAD
import timeit
#%%
Nx = 210
Nz = 68  
offset_max = 4000 
tEnd = 2.5 
dt = 0.004
t = np.arange(0,tEnd,dt)
frec = 3
Nt = tEnd/dt 
Nt= int(Nt)
dz = 25
dx = 25
vp_ori = np.zeros((Nx,Nz))
vp_ite = np.zeros((Nx,Nz))
Sx1 = 25
Sx2 = 185
Sx3 = 105
Sx4 = 65
Sx5 = 145
Sz = 3
Sz1 = 3
Sz2 = 3
Sz3 = 3
Sz4 = 3
Sz5 = 3
beta1 = 1
beta2 = 1
beta3 = 1
beta4 = 1
beta5 = 1
f1 = 0
f2 = 0
f3 = 0
f4 = 0
f5 = 0
grad1 = np.zeros((Nx*Nz,1))
grad2 = np.zeros((Nx*Nz,1)) 
grad3 = np.zeros((Nx*Nz,1)) 
grad4 = np.zeros((Nx*Nz,1)) 
grad5 = np.zeros((Nx*Nz,1)) 

a = (np.pi*frec)**2
t0 = 1
g1 = (-2*a*(t-t0)*np.exp(-a*(t-t0)**2)).T
g1 = np.reshape(g1,(np.size(g1),1))

# Creación del modelo de velocidades original

for iz in range(0,Nz):
    vp_ori[:,iz] = 2000+0.7*dz*(iz)

v11,v22 = np.shape(vp_ori)
X,Y = np.meshgrid(np.arange(0,v22,1), np.arange(0,v11,1))
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
# surf = ax.plot_surface(X,Y,vp_ori,cmap=cm.coolwarm,linewidth=5) 
# plt.show()

#%% Adquisición sobre el modelo original
start=timeit.default_timer()
Pt_obs1 = propagator(vp_ori, g1, Sx1, Sz, dx, dz, dt, 4000, frec)
Pt_obs2 = propagator(vp_ori, g1, Sx2, Sz, dx, dz, dt, 4000, frec)
Pt_obs3 = propagator(vp_ori, g1, Sx3, Sz, dx, dz, dt, 2125, frec)
Pt_obs4 = propagator(vp_ori, g1, Sx4, Sz, dx, dz, dt, 3125, frec)
Pt_obs5 = propagator(vp_ori, g1, Sx5, Sz, dx, dz, dt, 3125, frec)
stop=timeit.default_timer()
print('time:', stop-start) 
#%% Creación del modelo de velocidad inicial

for iz in range(0,Nz):
    vp_ite[:,Nz-iz-1] = 3172.5-0.6*dz*(iz)
    
v33,v44 = np.shape(vp_ite)
X1,Y1 = np.meshgrid(np.arange(0,v44,1), np.arange(0,v33,1))
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
# surf = ax.plot_surface(X1,Y1,vp_ite,cmap=cm.coolwarm,linewidth=5) 
X,Y = np.meshgrid(np.arange(0,v22,1), np.arange(0,v11,1))
# surf = ax.plot_surface(X,Y,vp_ori,cmap=cm.coolwarm,linewidth=5) 
# plt.show()

#%%
start=timeit.default_timer()
FF = 40
f = np.zeros(FF)
vp_ite1 = np.reshape(vp_ite,(Nx*Nz), order='F')
i=0
for i in range(0,FF):
# intentos=0
#while intentos<10:  
    #vp_ite1 = np.reshape(vp_ite,(Nx*Nz), order='F')
    # Cálculo del gradiente y función objetivo para la adquisición
    #1
    f1,grad1 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx1, Sz1, dx, dz, dt, 4000, frec, Pt_obs1[0], 1) # offset_max = 10000 k=6
    #2
    f2,grad2 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx2, Sz2, dx, dz, dt, 4000, frec, Pt_obs2[0], 2) # offset_max = 10000 k=6
    ##3
    f3,grad3 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx3, Sz3, dx, dz, dt, 2125, frec, Pt_obs3[0], 3) # offset_max = 10000 k=6
    #4
    f4,grad4 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx4, Sz4, dx, dz, dt, 3125, frec, Pt_obs4[0], 4) # offset_max = 10000 k=6
    #5
    f5,grad5 = FWI_GRAD(vp_ite1, Nx, Nz, Nt, g1, Sx5, Sz5, dx, dz, dt, 3125, frec, Pt_obs5[0], 5) # offset_max = 10000 k=6
    
    # función objetivo total y gradiente total
    ft = f1 + f2 + f3 + f4 + f5
            
    if i>1:
        if ft>f[i-1]:
            beta1= beta1*0.5
            beta2= beta2*0.5
            beta3= beta3*0.5
            beta4= beta4*0.5
            beta5= beta5*0.5
            f[i]=f[i-1]
            intentos= intentos+1
            print('valor rechazado=%f' % ft)
        else:
            f[i]=ft
            vp_ite=vp_ite1
            intentos=0
            print(f[i])
    else:
        f[i]=ft
        print(f[i])

    print("Iteración " + str(i))  
    #print('iteration=%d beta=%f' % i, beta1)
    # Aplicando el factor de escala
    grad = beta1*grad1 + beta2*grad2 + beta3*grad3 + beta4*grad4 + beta5*grad5
    
    
    # Actualizando el modelo de velocidades
    vp_ite1 = vp_ite + np.reshape(grad[:],(Nx,Nz),order='F')
    i=i+1
    #v55,v66 = np.shape(vp_ite)
    #X2,Y2 = np.meshgrid(np.arange(0,v66,1), np.arange(0,v55,1))
    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
    #surf = ax.plot_surface(X2,Y2,vp_ite,cmap=cm.nipy_spectral,linewidth=5) 
    #surf = ax.plot_surface(X,Y,vp_ori,cmap=cm.coolwarm,linewidth=5,alpha=0.8) 
    #plt.show()
stop=timeit.default_timer()
print('time:', stop-start)   
#%% Grafica del modelo de velocidades actualizado vs original
v55,v66 = np.shape(vp_ite)
X2,Y2 = np.meshgrid(np.arange(0,v66,1), np.arange(0,v55,1))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
surf = ax.plot_surface(X2,Y2,vp_ite,cmap=cm.nipy_spectral,linewidth=5) 
surf = ax.plot_surface(X,Y,vp_ori,cmap=cm.coolwarm,linewidth=5,alpha=0.8) 
plt.show()      