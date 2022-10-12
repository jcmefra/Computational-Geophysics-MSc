# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 17:01:29 2022

@author: juank
"""

# Very Fast simulated anhealing
# Implemented by Sergio Abreo
# Computational Geophysics
# October 2018

#from math import sqrt, 
from heuristics import punto2
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import statistics as st

intervalo1a=0.001
intervalo1b=0.02  
paso1=0.001      
intervalo2a=1
intervalo2b=10
paso2=1

N=2 # numero de parametros
M=1 # numero de observaciones

X,Y = np.meshgrid(np.arange(intervalo1a,intervalo1b,paso1), np.arange(intervalo2a,intervalo2b,paso2))
k1,k2=X.shape
Z=np.zeros((k1,k2))

for i in range(0, k1):
    for j in range(0, k2):
        a= np.array([X[i, j], Y[i,j]])
        Z[i][j]= punto2(a)

maxi=[intervalo1b, intervalo2b]
mini=[intervalo1a, intervalo2a]
T=np.array([10,8])
c=[0.05,0.05]

fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))

experimentos=100
resultados=np.array(np.zeros((experimentos,2)))
for l in range(experimentos):
    k=0
    m0=np.array([np.random.uniform(0.001,0.02),np.random.uniform(1,10)]) # Creación de temperaturas iniciales aleatorias para ambos parámetros
    Em0=punto2(m0)
    #plt.plot(m0[0],m0[1],Em0,'ro',label='Punto inicial')
    m1=m0.copy()
    while (T.min() > 0):
        k=k+1
        for i in range(100): # loop over a number of random moves/temperature
            for j in range(N*M):
                u=np.random.rand(1)
                yi=np.sign(u-0.5)*T[j]*((1+T[j])**(np.abs(2*u-1))-1)
                m1[j]=m1[j]+yi*(maxi[j]-mini[j])
                if m1[j]< mini[j]:
                    m1[j]= mini[j]
                if m1[j]>= maxi[j]:
                    m1[j]= maxi[j]

            Em1=punto2(m1)
            Delta_e=Em1-Em0
            P=np.exp(-Delta_e/T.min())
            if (Delta_e<=0):
                m0=m1.copy()
                Em0=Em1.copy()
                ax.plot(m0[0],m0[1],Em0,'b*')
            else:
                r=np.random.rand(1)
                if (P > r):
                    m0=m1.copy()
                    Em0=Em1.copy()
                    ax.plot(m0[0],m0[1],Em0,'b*')   
        for j in range(N*M):
            T[j]=T[j]*np.exp(-c[j]*k**(1/(M*N)))  
    resultados[l]=m0

prom_resultados=np.mean(resultados,0)
r_alpha=resultados[:,0]
r_n=resultados[:,1]
prom_alpha=np.mean(r_alpha)
prom_n=np.mean(r_n) 
mediana_alpha=st.median(r_alpha)
mediana_n=st.median(r_n)
varianza_alpha=np.var(r_alpha)
varianza_n=np.var(r_n)                      
        
surf = ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,alpha=0.3,linewidth=5)
ax.plot(m0[0],m0[1],Em0,'bo',label='Punto final')
plt.plot(m0[0],m0[1],Em0,'mo', label='Puntos iniciales') 
plt.show()        
#print(m0)
print(resultados)