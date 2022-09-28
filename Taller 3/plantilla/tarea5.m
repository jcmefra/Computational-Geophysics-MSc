% Tarea 5. Capítulo 3, segunda edición de Aster. Ejercicio 4.
% Plantilla elaborada por: Sergio A. Abreo Carrillo
% Curso de geofísica computaciona;

% Limpia las variables del espacio de trabajo
clear
clc
close all

%---------------------------------------
% Lectura de los dos archivos rowscan, colscan, diag1scan, diag2scan, std
%---------------------------------------

% Estructura que contiene toda la información del escaneo por filas de oeste a
% este

[x1 y1 x2 y2 t]=textread('rowscan','%f %f %f %f %f');
rowscan.x1=x1;
rowscan.y1=y1;
rowscan.x2=x2;
rowscan.y2=y2;
rowscan.t=t;
d1=rowscan.t;

% Estructura que contiene toda la información del escaneo por columnas de
% norte a sur

[x1 y1 x2 y2 t]=textread('colscan','%f %f %f %f %f');
colscan.x1=x1;
colscan.y1=y1;
colscan.x2=x2;
colscan.y2=y2;
colscan.t=t;
d2=colscan.t;

% Estructura que contiene toda la información del escaneo en diagonal
% west to east south to north

[x1 y1 x2 y2 t]=textread('diag1scan','%f %f %f %f %f');
diag1scan.x1=x1;
diag1scan.y1=y1;
diag1scan.x2=x2;
diag1scan.y2=y2;
diag1scan.t=t;
d3=diag1scan.t;

% Estructura que contiene toda la información del escaneo en diagonal
% west to east from north to south

[x1 y1 x2 y2 t]=textread('diag2scan','%f %f %f %f %f');
diag2scan.x1=x1;
diag2scan.y1=y1;
diag2scan.x2=x2;
diag2scan.y2=y2;
diag2scan.t=t;
d4=diag2scan.t;

% Desviación estándar de los datos
[sigma]=textread('std','%f');

% Para el primer punto el for va hasta 32
% para el segundo punto el for va hasta 94

for i=1:32 %32 94
    v(i)=1/sigma;
end

% Matriz de ponderaciones
W=diag(v,0);

% Primer punto
d=[d1; d2];

% Segundo punto
%d=[d1;d2;d3;d4];

dw=W*d;

%% Construcción de la matriz del sistema
k=sqrt(2);

% Escaneo por filas. West to east
G1=[ ones(1,16), zeros(1,240); 
   zeros(1,16), ones(1,16), zeros(1,224); 
   zeros(1,32), ones(1,16), zeros(1,208); 
   zeros(1,48), ones(1,16), zeros(1,192); 
   zeros(1,64), ones(1,16), zeros(1,176); 
   zeros(1,80), ones(1,16), zeros(1,160); 
   zeros(1,96), ones(1,16), zeros(1,144); 
  zeros(1,112), ones(1,16), zeros(1,128); 
  zeros(1,128), ones(1,16), zeros(1,112); 
  zeros(1,144), ones(1,16), zeros(1,96); 
  zeros(1,160), ones(1,16), zeros(1,80); 
  zeros(1,176), ones(1,16), zeros(1,64); 
  zeros(1,192), ones(1,16), zeros(1,48); 
  zeros(1,208), ones(1,16), zeros(1,32); 
  zeros(1,224), ones(1,16), zeros(1,16);
  zeros(1,240), ones(1,16)];

% Escaneo por columnas. North to south
G2=[1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15);
    0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14), 0, 1, zeros(1,14); 
    zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13), zeros(1,2), 1, zeros(1,13);
    zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12), zeros(1,3), 1, zeros(1,12);
    zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11), zeros(1,4), 1, zeros(1,11);
    zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10), zeros(1,5), 1, zeros(1,10);
    zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) , zeros(1,6), 1, zeros(1,9) ;
    zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) , zeros(1,7), 1, zeros(1,8) ;
    zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) , zeros(1,8), 1, zeros(1,7) ;
    zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) , zeros(1,9), 1, zeros(1,6) ;
    zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5), zeros(1,10), 1, zeros(1,5);
    zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4), zeros(1,11), 1, zeros(1,4);
    zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3), zeros(1,12), 1, zeros(1,3);
    zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2), zeros(1,13), 1, zeros(1,2);
    zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1), zeros(1,14), 1, zeros(1,1);
    zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1, zeros(1,15), 1];

% Escaneo en diagonal West to East, South to North
G3=[ k, zeros(1,255);
     0, k, zeros(1,14), k, zeros(1,239);
     0, 0, k, zeros(1,13), 0, k, zeros(1,14), k, zeros(1,223); 
     0, 0, 0, k, zeros(1,12), 0, 0, k, zeros(1,13), 0, k, zeros(1,14), k, zeros(1,207);
     zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,191);
     zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,175);
     zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,159);
     zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,143);
     zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,127);
     zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,111);
     zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,95);
     zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,79);
     zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,63);
     zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,47);
     zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,31);
     zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14), k, zeros(1,15);
     zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13), zeros(1,1), k, zeros(1,14);
     zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12), zeros(1,2), k, zeros(1,13);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11), zeros(1,3), k, zeros(1,12);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10), zeros(1,4), k, zeros(1,11);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9), zeros(1,5), k, zeros(1,10);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8), zeros(1,6), k, zeros(1,9);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7), zeros(1,7), k, zeros(1,8);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6), zeros(1,8), k, zeros(1,7);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5), zeros(1,9), k, zeros(1,6);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4), zeros(1,10), k, zeros(1,5);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3), zeros(1,11), k, zeros(1,4);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2), zeros(1,12), k, zeros(1,3);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1), zeros(1,13), k, zeros(1,2);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k, zeros(1,14), k, zeros(1,1);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,15), k];

% Escaneo en diagonal West to East, North to South
G4=[ zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5);
     zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4);
     zeros(1,16), zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3);
     zeros(1,16), zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2);
     zeros(1,16), k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1);
     k, zeros(1,15), zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k;
     zeros(1,1), k, zeros(1,14), zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16);
      zeros(1,2), k, zeros(1,13), zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16);
      zeros(1,3), k, zeros(1,12), zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,4), k, zeros(1,11), zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,5), k, zeros(1,10), zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,6), k, zeros(1,9), zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,7), k, zeros(1,8), zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,8), k, zeros(1,7), zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,9), k, zeros(1,6), zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,10), k, zeros(1,5), zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,11), k, zeros(1,4), zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,12), k, zeros(1,3), zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,13), k, zeros(1,2), zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,14), k, zeros(1,1), zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16);
      zeros(1,15), k, zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16), zeros(1,16) ];

% primer punto
    %Matriz G para el primer punto tomando filas y columnas
    G12 = [G1;G2];
    rank_1 = rank(G12);
    
    %Información de USV
    [U, S, V] = svd(G12);
    %USVp
    [Up, Sp, Vp] = svdsketch(G12,0.01);

    %G Daga
    G_dagger = Vp*inv(Sp)*Up';
    m_dagger = G_dagger*d;
    max_m = max(m_dagger)
    min_m = min(m_dagger)


    %Resolución 
    Rm = G_dagger * G12;
    Vector_diag = diag(Rm);
    Rmdiag = diag(Vector_diag,0);


%% Visualizando 

% el bloque
figure
imagesc(reshape(m_dagger,16,16))
colormap(gray)
title('Inversión del bloque')
colorbar

figure
contour(reshape(m_dagger,16,16))
set(gca, 'YDir','reverse')
title('Inversión del bloque')
colorbar

% la matriz de resolución
figure
mesh(Rm)
title('Matriz de resolución')
colorbar

figure
imagesc(Rmdiag)
colormap(gray)
title('Elementos de la diagonal de la matriz de resolución')
colorbar

figure
mesh(Rmdiag)
title('Elementos de la diagonal de la matriz de resolución')
colorbar

%Valores singulares
Diag_S = diag(S);

%Espacio Nulo

Vo=V(:,99);
Vo_Matrix = reshape(Vo,16,16);
Vo_diag = diag(Vo_Matrix);
sum_Vo = sum(Vo_Matrix,'all')
sumdiag = sum(Vo_diag,"all")
figure
imagesc(reshape(Vo,16,16))
colormap(gray)
colorbar
