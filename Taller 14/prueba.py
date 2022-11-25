from scipy import optimize
import numpy as np
from funcfwi import funcfwi
from gradfwi import gradfwi
import matplotlib.pyplot as plt
from matplotlib import cm
import timeit

# funcion base cuadratica
#def funcion(x):
#    y = x**2
#    return y

# derivada de la funcion base
#def fprima(x):
#    y = 2*x
#    return y

#x0=5
#minimo=optimize.fmin_l_bfgs_b(funcion,x0,fprime=fprima) 

Nx = 210
Nz = 68
dz = 25
vp_ite = np.zeros((Nx,Nz))

# Modelo inicial
for iz in range(0,Nz):
    vp_ite[:,Nz-iz-1] = 3172.5-0.6*dz*(iz)

vp_ite1 = np.reshape(vp_ite,(Nx*Nz), order='F')

v11,v22 = np.shape(vp_ite)
X,Y = np.meshgrid(np.arange(0,v22,1), np.arange(0,v11,1))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
surf = ax.plot_surface(X,Y,vp_ite,cmap=cm.coolwarm,linewidth=5) 
plt.show()

#inicial=funcfwi(vp_ite)
#print(inicial)

#ginicial=gradfwi(vp_ite)
#print(ginicial)

#gview=np.reshape(ginicial[:],(Nx,Nz),order='F')
#v33,v44 = np.shape(gview)
#X1,Y1 = np.meshgrid(np.arange(0,v44,1), np.arange(0,v33,1))
#fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
#surf = ax.plot_surface(X1,Y1,gview,cmap=cm.coolwarm,linewidth=5) 
#plt.show()

start = timeit.default_timer()
minimo=optimize.fmin_l_bfgs_b(funcfwi,vp_ite1,fprime=gradfwi,iprint=1) 
#minimo = optimize.minimize(funcfwi, vp_ite1, method='L-BFGS-B', jac=gradfwi)
stop = timeit.default_timer()

print('Tiempo L-BFGS-B:',stop-start)
print(minimo)

#mview=np.reshape(minimo[0],(Nx,Nz),order='F')
mview=np.reshape(minimo.x,(Nx,Nz),order='F')
v33,v44 = np.shape(mview)
X1,Y1 = np.meshgrid(np.arange(0,v44,1), np.arange(0,v33,1))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
surf = ax.plot_surface(X1,Y1,mview,cmap=cm.nipy_spectral,linewidth=5)  
#surf = ax.plot_surface(X,Y,vp_ori,cmap=cm.coolwarm,linewidth=5,alpha=0.8)
plt.show()