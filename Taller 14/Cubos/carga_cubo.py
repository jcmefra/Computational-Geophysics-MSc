import numpy as np
import matplotlib.pyplot as plt
from drawnow import drawnow

#Set the dimensions
Nx=200;
Nz=200;


# Adjust the path to your computer
fid = open ('Cubos/h_vel.bin','rb'); L=50;
#fid = open ('Cubos/h_rd.bin','rb'); L=50;
#fid = open ('Cubos/gk_c.bin','rb'); L=10;
#fid = open ('Cubos/mk_c.bin','rb'); L=10;
prueba = np.fromfile (fid, dtype = np.float32);
fid.close();
prueba = prueba.reshape(int((len(prueba))/Nz), Nz);
prueba = np.transpose(prueba);

video = np.zeros([Nz, Nx, L])
for k in range (0, L):
    for j in range (0, Nx):
        for i in range (0, Nz):
                video [i, j, k] = (prueba[i,j+(k-1)*Nx]); 


def make_fig():
    plt.imshow(video[:,:,i])
    plt.title('Iteration ' + str(i+1))
    plt.ylabel('Depth (m)')
    plt.xlabel('Distance (m)')

for i in range(0, L):
    drawnow(make_fig)
plt.show()

# for i=1:L
#     imagesc(video(:,:,i)),colorbar
#     xlabel('Distance (m)');
#     ylabel('Depth (m)');
#     str=sprintf('Iteracion %d',i);title(str);
#     pause(0.25)
# end