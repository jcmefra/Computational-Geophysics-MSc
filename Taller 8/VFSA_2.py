from tqdm import tqdm
from heuristics import punto2
import numpy as np
import matplotlib.pyplot as plt

alpha_0 = []
n_0 = []
muestreo_plot = []
intervalo1a=0.001
intervalo1b=0.02  
paso1=0.001       
intervalo2a=1
intervalo2b=10
paso2=0.1

maxi=[intervalo1b, intervalo2b]
mini=[intervalo1a, intervalo2a]
T=np.array([10 ,8])
c=[0.005,0.05]

N=2 # numero de parametros
M=1 # numero de observaciones
ensayos = 100
for ensayos in tqdm(range(100)):
    m0=np.array([np.random.uniform(0.001,0.02),np.random.uniform(1,10)])
    m0=m0.ravel()
    Em0=punto2(m0)
    m1=[m0[0],m0[1]]
    k=0
    while (T[0] > 0 and T[1] > 0): # loop over the temperature T
        k+=1
        for i in range (100): # loop over a number of random moves/temperature
            for j in range(N*M):
                u=np.random.rand(1)
                yi=np.sign(u-0.5)*T[j]*((1+T[j])**(np.abs(2*u-1))-1)
                m1[j]=m1[j]+yi*(maxi[j]-mini[j])
                if (m1[j]< mini[j]):
                    m1[j]= mini[j]
                elif (m1[j]>= maxi[j]):
                    m1[j]= maxi[j]
            Em1=punto2(m1)
            Delta_e=Em1-Em0
            P=np.exp(-Delta_e/T[0])
            if (Delta_e<=0):
                m0=m1.copy()
                Em0=Em1.copy()
            else:
                r=np.random.rand(1)
                if (P > r):
                    m0=m1.copy()
                    Em0=Em1.copy()
            muestreo_plot.append([m0[0],m0[1],Em0])
        for j in range(N*M):
            T[j]=T[j]*np.exp(-c[j]*k**(1/(M*N)))  
    alpha_0.append(m0[0])
    n_0.append(m0[1])         

print ("alpha_0 es igual a", np.mean(alpha_0))
print ("n_0 es igual a", np.mean(n_0))
print ("Varianza alpha", np.var(alpha_0))
print ("Varianza n", np.var(n_0))
print ("Mediana alpha", np.median(alpha_0))
print ("Mediana n", np.median(n_0))

print (len(muestreo_plot))
muestreo_plot_arr = np.array(muestreo_plot)
fig,ax = plt.subplots(subplot_kw = {"projection":"3d"}, figsize = (8,8))
ax.scatter(muestreo_plot_arr[:,0], muestreo_plot_arr[:,1], muestreo_plot_arr[:,2], alpha =0.5)
ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
ax.legend(loc='upper left')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()     

plt.title('Histograma Alpha')
plt.hist(alpha_0)
plt.show()

plt.title('Histograma n')
plt.hist(n_0)
plt.show()