# Very Fast simulated anhealing
# Implemented by Sergio Abreo
# Computational Geophysics
# October 2018

from heuristics import boha1
from heuristics import drop
from heuristics import matyas
from heuristics import ackley 
from matplotlib import cm
import numpy as np 
import matplotlib.pyplot as plt

N=2 # numero de parametros
M=1 # numero de observaciones

intervalo=10  #32 ackl, 100 boha1,   5 drop, 10 matya
paso=1      # 1 ackl,   1 boha1, 0.5 drop,  1 matya

X,Y = np.meshgrid(np.arange(-intervalo,intervalo+1,paso), np.arange(-intervalo,intervalo+1,paso))

k=np.shape(X)
Z=np.empty(np.shape(X))
k=k[0]

for i in range(k):
    for j in range(k):
        #Z[i,j]=ackley( [ X[i,j] , Y[i,j] ] )
        #Z[i,j]=boha1( [ X[i,j] , Y[i,j] ] )
        #Z[i,j]=drop( [ X[i,j] , Y[i,j] ] )
        Z[i,j]=matyas( [ X[i,j] , Y[i,j] ] )
        
T=np.zeros((1,2))
c=np.zeros((1,2))
maxi=np.zeros((1,2))
mini=np.zeros((1,2))
m0=np.zeros((1,2))
T=T.flatten().tolist() 
c=c.flatten().tolist()    
maxi=maxi.flatten().tolist()
mini=mini.flatten().tolist()   
m0=m0.flatten().tolist() 

To=[15,15]
Co=[0.05,0.05]
for j in range((M*N)):
    T[j]=To[j] # Temperatura inicial
    c[j]= Co[j] # coeficientes de cada parametro
    maxi[j]=intervalo # maximo valor del parametro
    mini[j]=-intervalo # minimo valor del parametro
    m0[j]=float(np.random.rand(1)*2 -1)*intervalo # 

#Em0=ackley(m0)
#Em0=boha1(m0)
#Em0=drop(m0)
Em0=matyas(m0)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
ax.plot(m0[0],m0[1],Em0,'ro',label='Punto inicial')
k=0
m1=[m0[0],m0[1]]
itt=0
while (T[0] > 0 and T[1] > 0):
    k=k+1
    for i in range(100): # loop over a number of random moves/temperature
        for j in range(N*M):
            u=float(np.random.rand(1))
            yi=float(np.sign(u-0.5)*T[j]*((1+T[j])**(np.abs(2*u-1))-1))
            m1[j]=float(m1[j]+yi*(maxi[j]-mini[j]))
            if (m1[j]< mini[j]):
                m1[j]= mini[j]
            elif (m1[j]>= maxi[j]):
                m1[j]= maxi[j]

        #Em1=ackley(m1)
        #Em1=boha1(m1)
        #Em1=drop(m1)
        Em1=matyas(m1)
        Delta_e=Em1-Em0
        P=np.exp(-Delta_e/T[0])
        if (Delta_e<=0):
            m0=[m1[0],m1[1]]
            Em0=Em1
            ax.plot(m0[0],m0[1],Em0,'b*')
            itt=itt+1
        else:
            r=float(np.random.rand(1))
            if (P > r):
                m0=[m1[0],m1[1]]
                Em0=Em1
                ax.plot(m0[0],m0[1],Em0,'b*') 
                itt=itt+1
    for j in range(N*M):
        T[j]=T[j]*np.exp(-c[j]*k**(1/(M*N)))
               
surf = ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,alpha=0.3,linewidth=5)
ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
plt.plot(m0[0],m0[1],Em0,'mo') 
plt.show()          
print(m0)