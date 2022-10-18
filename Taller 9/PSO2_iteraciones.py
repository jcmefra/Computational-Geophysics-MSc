from heuristics import punto2
import numpy as np
import numba as nb
from tqdm import tqdm
from matplotlib import cm
import matplotlib.pyplot as plt

#start=timeit.default_timer()

alpha_0 = []
n_0 = []
muestreo_plot = []
intervalo1a=0.001
intervalo1b=0.02  
paso1=0.001       
intervalo2a=1
intervalo2b=10
paso2=0.1

Xf,Yf  = np.meshgrid(np.arange(intervalo1a,intervalo1b,paso1), np.arange(intervalo2a,intervalo2b,paso2))
k1,k2=Xf.shape
Zf=np.empty(np.shape(Xf))

for i in range(0, k1):
    for j in range(0, k2):
        a= np.array( [ Xf[i,j] , Yf[i,j] ] )
        Zf[i][j]= punto2(a)
        
particulas=20
iteraciones=100

V=np.zeros((particulas,2))

#parámetros PSO
w=0.729 
c1=2.05
c2=2.05
X=np.zeros((particulas,2))
Z=np.zeros((1,particulas))

experimentos=50
resultados=np.array(np.zeros((experimentos,2)))
a = 0
fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
for l in tqdm(range(experimentos)):
    a+=1
    k=0
    for i in range(particulas):
        X[i,:]=np.array([np.random.uniform(intervalo1a,intervalo1b),np.random.uniform(intervalo2a,intervalo2b)])
        Z[0,i]=punto2(X[i,:])

    # 2) Encontrar valores de Pi y Pg
    Pi=X
    Y=Z.min()
    I=Z.argmin()
    Pg=X[I,:]

    #ax.scatter(X[:,0], X[:,1], Z, s = 60, color = 'b',label='Enjambre inicial')  
    if a == 1:
        surf = ax.plot_surface(Xf,Yf,Zf,cmap=cm.coolwarm,alpha=0.3,linewidth=5) 

    for j in range(iteraciones):
        #3) Actualizar posicion y velocidad
        for i in range(particulas):
            r1=np.random.rand(1)
            r2=np.random.rand(1)
            V[i,:]=w*V[i,:] + c1*r1*(Pi[i,:]-X[i,:]) + c2*r2*(Pg-X[i,:])
            Xt=X[i,:] + V[i,:]
            if Xt[0]<intervalo1a:
                Xt[0]=intervalo1a
            if Xt[0]>intervalo1b:
                Xt[0]=intervalo1b
            if Xt[1]<intervalo2a:
                Xt[1]=intervalo2a
            if Xt[1]>intervalo2b:
                Xt[1]=intervalo2b
            
            Yt=punto2(Xt)
                    
            # Evalua si el nuevo valor es mejor para actualizar la posicion de la particula
            if Z[0,i] > Yt:
                X[i,:] = Xt
            
            # Obtener nuevo Pi
            Zp=punto2(X[i,:])
            
            if Z[0,i]> Zp:
                Pi[i,:]=X[i,:]
                Z[0,i]=Zp
            else:
                Pi[i]=Pi[i]
                
            # Obtener el nuevo Pg
            Yg=punto2(Pg)

            if Yg > Zp:
                Pg=X[i,:]
            else:
                Pg=Pg
    ax.scatter(X[:,0], X[:,1], Z, s = 60, alpha =0.25, color = 'b',label='Enjambre final') 
    ax.set_title('Experimento={}'.format(l+1))
    ax.set_xlim([intervalo1a,intervalo1b])
    ax.set_ylim([intervalo2a,intervalo2b])
    ax.set_xlabel('α')
    ax.set_ylabel('n')
    ax.set_zlabel('Z')
    if a == 1:
        ax.legend(loc='upper left')
    fig.canvas.draw()
    fig.canvas.flush_events()       

    resultados[l]=Pg
    alpha_0.append(resultados[:,0])
    n_0.append(resultados[:,1])

#surf = ax.plot_surface(Xf,Yf,Zf,cmap=cm.coolwarm,alpha=0.3,linewidth=5) 

#stop=timeit.default_timer()
#print('Time: ', stop - start)
# print(Pg.flatten().tolist())
print (resultados)
#%% Estadísticas
print ("alpha_0 es igual a", np.mean(alpha_0))
print ("n_0 es igual a", np.mean(n_0))
print ("Varianza alpha", np.var(alpha_0))
print ("Varianza n", np.var(n_0))
print ("Mediana alpha", np.median(alpha_0))
print ("Mediana n", np.median(n_0))

# plt.title('Histograma Alpha')
# plt.hist(alpha_0)
# plt.show()

# plt.title('Histograma n')
# plt.hist(n_0)
# plt.show()

# print (len(muestreo_plot))
# muestreo_plot_arr = np.array(muestreo_plot)
# fig,ax = plt.subplots(subplot_kw = {"projection":"3d"}, figsize = (8,8))
# ax.scatter(muestreo_plot_arr[:,0], muestreo_plot_arr[:,1], muestreo_plot_arr[:,2], alpha =0.5)
# ax.plot(X[:,0], X[:,1], Z, s = 60, color = 'r',label='Enjambre final')
# ax.legend(loc='upper left')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# plt.show()  

# prom_resultados=np.mean(resultados,0)

