%Computational Geophysics Course
%Dr.Ing. Sergio Abreo.
%Loading the cubes of optimal steps
%November, 2022

clear all
clc 
close all

%Set the dimensions
Nx=200;
Nz=200;

% Adjust the path to your computer
%fid = fopen('h_vel.bin','rb'); L=50;
fid = fopen('h_rd.bin','rb'); L=50;
%fid = fopen('gk_c.bin','rb'); L=10;
%fid = fopen('mk_c.bin','rb'); L=10;
prueba = fread(fid,'float32');
fclose(fid);
[prueba,pad_f1] = vec2mat(prueba,Nz);
prueba = prueba';

for k=1:L
    for j=1:Nx
        for i=1:Nz
            video(i,j,k)=prueba(i,j+(k-1)*Nx); 
        end
    end
end
%
figure(1)
for i=1:L
    imagesc(video(:,:,i)),colorbar
    xlabel('Distance (m)');
    ylabel('Depth (m)');
    str=sprintf('Iteracion %d',i);title(str);
    pause(0.25)
end