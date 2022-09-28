import heuristics
from heuristics import boha1
from heuristics import drop
from heuristics import matyas
from heuristics import ackley 
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


xx=[-32, 32]
xx=np.array(xx)
matyas(xx)
ackley(xx)
boha1(xx)
drop(xx)
intervalo=32;  # 32 ackl, 100 boha1,   5 drop, 10 matyas
paso=1;        # 1 ackl,   1 boha1, 0.5 drop,  1 matyas

X,Y = np.meshgrid(np.arange(-intervalo,intervalo+1,paso), np.arange(-intervalo,intervalo+1,paso))
k=len(X)
Z=np.zeros((k,k))

for i in range(0, k):
    for j in range(0, k):
        a= np.array([X[i, j], Y[i,j]])
        #Z[i][j]= matyas(a)
        Z[i][j]= ackley(a)
        #Z[i][j]= drop(a)
        #Z[i][j]= boha1(a)

T=10; 
m0=(np.random.rand(1,2)*2 -1)*intervalo
m0=m0.ravel()
Em0=ackley(m0)
#Em0=boha1(m0)
#Em0=drop(m0)
#Em0=matyas(m0)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
#view(-17,18)#ackley
#view(-69,67)#boha1
#view(-161,29)#drop
#ax.view_init(30,angle)#matya

plt.plot(m0[0],m0[1],Em0,'ro',label='Punto inicial')

while T > 0: # loop over the temperature T
    for i in range (100): # loop over a number of random moves/temperature
        m1=(np.random.rand(1,2)*2 -1)*intervalo
        m1=m1.ravel()
        Em1=ackley(m1)
        #Em1=boha1(m1)
        #Em1=drop(m1)
        #Em1=matyas(m1)
        Delta_e=Em1-Em0
        P=np.exp(-Delta_e/T)
        if (Delta_e<=0):
           m0=m1
           Em0=Em1
           plt.plot(m0[0],m0[1],Em0,'b*')
        else:
            r=np.random.rand(1)
            if (P > r):
                m0=m1
                Em0=Em1
                plt.plot(m0[0],m0[1],Em0,'b*') 
    if(T >= 2):
        T=T-1
    else:
        T=T-0.2

print(m0)       
surf = ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,alpha=0.3,linewidth=5)
ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
plt.plot(m0[0],m0[1],Em0,'mo') 
plt.show()        