#import heuristics
from tqdm import tqdm
from heuristics import punto2
import numpy as np
import matplotlib.pyplot as plt

alpha_0 = []
n_0 = []
muestreo_plot = []
w = 0
for w in tqdm(range(0,100)):
    w+=1
    intervalo1a=0.001
    intervalo1b=0.02  
    paso1=0.001       
    intervalo2a=1
    intervalo2b=10
    paso2=1

    T=10; 
    temp1=(np.random.rand(1,1)*0.019 +0.001) #[0.01-0.02]
    temp2=(np.random.rand(1,1)*9 +1) #[1-10] 
    m0=np.array([temp1, temp2])
    m0=m0.ravel()
    Em0=punto2(m0)

    while T > 0: # loop over the temperature T
        for i in range (100): # loop over a number of random moves/temperature
            temp1=(np.random.rand(1,1)*0.019 +0.001) 
            temp2=(np.random.rand(1,1)*9 +1)
            m1=np.array([temp1, temp2])
            m1=m1.ravel()
            Em1=punto2(m1)
            Delta_e=Em1-Em0
            P=np.exp(-Delta_e/T)
            if (Delta_e<=0):
                m0=m1
                Em0=Em1
            else:
                r=np.random.rand(1)
                if (P > r):
                    m0=m1
                    Em0=Em1
            muestreo_plot.append([m0[0],m0[1],Em0])
        if(T >= 2):
            T=T-1
        else:
            T=T-0.2

    alpha_0.append(m0[0])
    n_0.append(m0[1])         

print ("alpha_0 es igual a", np.mean(alpha_0))
print ("n_0 es igual a", np.mean(n_0))
print ("Varianza alpha", np.var(alpha_0))
print ("varianza n", np.var(n_0))
print (len(muestreo_plot))
muestreo_plot_arr = np.array(muestreo_plot)
fig,ax = plt.subplots(subplot_kw = {"projection":"3d"}, figsize = (8,8))
ax.plot_trisurf(muestreo_plot_arr[:,0], muestreo_plot_arr[:,1], muestreo_plot_arr[:,2], edgecolor='none', alpha =.0025, antialiased = True)
ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
ax.legend(loc='upper left')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()        
