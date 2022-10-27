import numpy as np
from math import pi
from math import exp
import matplotlib.pyplot as plt
import timeit
from numba import jit
import ricker as rck
#%%

Tout=200 #Propagation time (rebotes fronteras no naturales)
Nx=200   #Number of point in the x direction
Nz=200   #Number of point in the z direction
c=1200.0 #Background velocity.
dh=0.75   #spatial step (resolucion espacial)
#dt = 0.0011
dt=dh/(c*np.sqrt(2))
G=c*dt/dh #Courant parameter.
fq=100 #Play for dispersion numerica from 10(ondicula pequeña) to 100(ondicula aumenta)
Sx=100
Sz=100

fuente=np.zeros(Tout)
for tstep in range(Tout):
    fuente[tstep]=rck.ricker(tstep, dt, fq) #source

@jit
def Propaga(Tout, Nx, Nz, dt, fq, Sx, Sz, fuente):
    P1 = np.zeros((Nz,Nx)) #pasado
    P2 = np.zeros((Nz,Nx)) #presente
    P3 = np.zeros((Nz,Nx)) #futuro 

    video = np.zeros((Nz,Nx,Tout))    

    # for tstep in range(Tout):
    #     fuente[tstep]=ricker(tstep, dt, fq) #source

    
    for t in range(Tout):
        for a in range(1,Nz-1):
            for b in range(1,Nx-1):
                P3[a,b] = (2-4*G*G)*P2[a,b]+G*G*(P2[a+1,b]+P2[a-1,b]+P2[a,b+1]+P2[a,b-1])-P1[a,b] 

        for c in range(1,Nz-1):
            for d in range(1,Nx-1):
                P1[c,d] = P2[c,d]
                P2[c,d] = P3[c,d]
                
        P2[Sz,Sx]=P2[Sz,Sx]+fuente[t] #inject the source
            
        for e in range(Nz):
            for f in range(Nx):
                video[e,f,t] = P3[e,f]
    
    return video

start=timeit.default_timer()
video=Propaga(Tout,Nx,Nz,dt,fq,Sx,Sz,fuente)
stop=timeit.default_timer()
print('Time: ', stop - start)

#Gráfica 1

paso = 1
Pasoo = np.arange(0,Tout,paso)
jj = np.shape(Pasoo)
fig, ax = plt.subplots(figsize=(10,10))
t=0
    
for t in range(jj[0]):
   plt.figure(1); plt.clf()
   plt.imshow(video[:,:,t])
   plt.title('Number ' + str(t+1))
   ax.set_ylabel('Depth (m)')
   plt.pause(0.1)

plt.show()
print(dt)