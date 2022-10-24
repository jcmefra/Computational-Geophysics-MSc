import numpy as np
import matplotlib.pyplot as plt
import timeit
from numba import jit

#%%

Tout=100 #Propagation time
Nx=200   #Number of point in the x direction
Nz=200   #Number of point in the z direction
c=1200.0 #Background velocity.
dh=2.0   #spatial step
dt=dh/(c*np.sqrt(2))
G=c*dt/dh #Courant parameter.
@jit
def Propaga(Tout, Nx, Nz):
    P1 = np.zeros((Nz,Nx))
    P2 = np.zeros((Nz,Nx))
    P3 = np.zeros((Nz,Nx))

    alpha = 0.15

    x0 = dh*Nx/2
    z0 = dh*Nz/3

    video = np.zeros((Nz,Nx,Tout))

    for i in range(Nz):
        for j in range(Nx):
            r2 = (dh*i-z0)*(dh*i-z0)+(dh*j-x0)*(dh*j-x0)
            P2[i,j] = np.sin(1-alpha*r2)*np.exp(-alpha*r2) #Source

    for t in range(Tout):
        for a in range(1,Nz-1):
            for b in range(1,Nx-1):
                P3[a,b] = (2-4*G*G)*P2[a,b]+G*G*(P2[a+1,b]+P2[a-1,b]+P2[a,b+1]+P2[a,b-1])-P1[a,b]

        for c in range(1,Nz-1):
            for d in range(1,Nx-1):
                P1[c,d] = P2[c,d]
                P2[c,d] = P3[c,d]
            
        for e in range(Nz):
            for f in range(Nx):
                video[e,f,t] = P3[e,f]
    
    return video

start=timeit.default_timer()
video=Propaga(100,200,200)
stop=timeit.default_timer()
print('Time: ', stop - start)

#Gráfica 1

paso = 1
Pasoo = np.arange(0,Tout,paso)
jj = np.shape(Pasoo)
fig, ax = plt.subplots(figsize=(10,10))
t=100
for t in range(jj[0]):
   plt.figure(1); plt.clf()
   plt.imshow(video[:,:,t])
   plt.title('Number ' + str(t+1))
   plt.pause(0.1)

plt.show()
