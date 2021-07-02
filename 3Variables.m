%% Carlos Chamizo Cano
%% indice
%{
% Introduccion								
% declaracion de la funcion
% grafico de la funcion 				
% curvas de nivel		
% grafica de la funcion junto con sus curvas de nivel
% grafica junto a curvas de nivel
% curva de nivel en una altura determinada
% valor de la funcion en un punto (x,y)
% limites
% Limites reiterados en un punt x_0 y_0
% continuidad de la funcion
% deribabilidad
% Derivadas parciales primeras 
% Representacion grafica de las primeras derivadas parciales 
% grafica curvas de nivel y gradiente
% Representacion del modulo del gradiente
% Grafica de la Hessiana por componentes
% campo de gradientes y de la hessiana
% Representacion de la norma del gradiente 
% Restricciones de la funcion sobre la frontera 
% Aproximacion de Taylor  
% Polinomio de Taylor de orden n centrado en [0,0]
% Apartado Errores asbolutos y relativos en el punto de la aproximacion
% optimizacion
% Puntos criticos
% Clasificacion mediante el criterio
% diferenciabilidad
% Aproximacion lineal en los puntos criticos
% superficie y de los planos tangentes
% aparatado creativo
%}
%% introduccion
%{
la funcion de esta practica es una "representacion" de la conservacion de
la energia en el universo, es decir si hay un porcentaje de materia,
aparece de misma cantidad esta de forma contrapuesta de antimateria. por lo es una 
una forma de representar que la energia ni se crea ni se destruye. Esta
practica esta menos trabajada que la anterior pero he tratado de dar la
mayor funcionalidad posible.
%}
%% declaracion de la funcion
clear  
syms x y z
%f(x,y)=1.3*exp(-(-cos(x-2.7)^2/0.32)+(sin(y+3)^2/0.32))*cos(x^2);
%la idea era usar esta funcion pero matlab no soportaba  

f(x,y)=y*1.2*exp(-(x^2+y^2))*((-x^3+y^2)^3)*(y * x)^2;  %funcion general usada
%f(x,y)=2*exp(-x^2-y^2)     %funcion con punto de silla para demostrar que
                            %se pueden detectar puntos de silla
f_n=matlabFunction(f,'Vars',{x , y});
%% grafico de la funcion
close all
fsurf(f)
%% curvas de nivel
%% grafica de la funcion junto con sus curvas de nivel
close all
figure
fcontour(f) %curvas de nivel (funciona tambien para z numerica)
%% grafica junto a curvas de nivel
close all
fsurf(f,[-5 5],'ShowContours','on')
%% curva de nivel en una altura determinada
C=1; %altura donde se desea encontrar la curva de nivel
curva=f-C; %resolvemos la ecuacion f=c cuya solucion es la curva deseada
subplot(2,1,1)
fimplicit(curva)
title('curvas de nivel en el z= 1')
subplot(2,1,2)
fcontour(f_n)
title('curvas de nivel de la superficie')
%% valor de la funcion en un punto (x,y)
x_0=0;
y_0=0;
punto=[x_0,y_0];
disp('la funcion vale z=');
disp(subs(f,[x,y],[0,0])) %valor de la funcion simbolica en le punto (0,2)
%% limites
%% Limites reiterados en un punt x_0 y_0
%{
se calcularan los distintos limites de la funcion en un punto que se decida
y despues se comprobará la existencia de los limites dobles
%}
x_0=0; y_0=0;
P=[x_0,y_0];
L1=limit(limit(f,x,x_0),y,y_0);  
L2=limit(limit(f,y,y_0),x,x_0);  
lim1=eval(L1);
lim2=eval(L2);
disp("/////limites reiterados///////")
disp(lim1)
disp(lim2)
%comprobacion de si existe el limite doble
if lim1==lim2
    disp("existe limite doble en z=");
    disp(lim1);
else
    disp("no existe limite doble ");
end
% Limites direccionales
dir_x=y^2;
dir_y=x^2;
y_dir=limit(subs(f,[x,y],[dir_x,y]),y,y_0);
x_dir=limit(subs(f,[x,y],[x,dir_y]),x,x_0);
lim1=eval(y_dir);
lim2=eval(x_dir);
disp("/////limites direccionales///////")
disp(lim1)
disp(lim2)
%comprobacion de si existe el limite doble
if lim1==lim2
    disp("existe limite doble en z=");
    disp(lim1);
else
    disp("no existe limite doble ");
end
% Limites radiales en (0,0)
syms k
y_r=y_0+k*(x-x_0);
x_r=x_0+k*(y-y_0);
LRX=limit(limit(f,y,y_r),x,x_0);
LRY=limit(limit(f,x,x_r),y,y_0);

disp("/////limites radiales///////")
disp(lim1)
disp(lim2)
%comprobacion de si existe el limite doble
if lim1==lim2
    disp("existe limite doble en z=");
    disp(lim1);
else
    disp("no existe limite doble ");
end
%% continuidad de la funcion
figure
fsurf(f,[-5 5 -5 5])
[num,den] = numden(f(x,y));
raices_deny = solve(den,y);
raices_denx = solve(den,x);
%{
la idea es comprobar que no haya raices del denominador de la funcion, y en
caso afirmativo significa que para ese valor no existe el dominio y por lo
tanto no es  continuo
%}
if isempty(raices_deny)==0 || isempty(raices_denx)==0
if isempty(raices_deny)==0
    disp("no es continua en y=")
    disp(raices_deny)
end
if isempty(raices_denx)==0
    disp("no es continua en x=")
    disp(raices_denx)
end
else 
disp("es continua en todo su dominio")
end
%% deribabilidad
%% Derivadas parciales primeras 
dfx=diff(f,x);  %derivada parcial de x
dfy=diff(f,y);  %derivada parcial de y
for i=1:2
    if i==1
        [num,den] = numden(dfx(x,y));
        raices_deny = solve(den,y);
        raices_denx = solve(den,x);
    else
        [num,den] = numden(dfy(x,y));
        raices_deny = solve(den,y);
        raices_denx = solve(den,x);
    end
%{
la idea es comprobar que no haya raices del denominador de la funcion, y en
caso afirmativo significa que para ese valor no existe el dominio y por lo
tanto no es  continuo
%}
if isempty(raices_deny)==0 || isempty(raices_denx)==0
if isempty(raices_deny)==0
    disp("no es deivable para los puntos con y =")
    disp(raices_deny)
end
if isempty(raices_denx)==0
    disp("no es derivable para los puntos con x =")
    disp(raices_denx)
end
else 
    if i==2
disp("es derivable en todo su dominio")
    end
end
end
%% Representacion grafica de las primeras derivadas parciales 
close all
figure
subplot(2,1,1)
fsurf(dfx,'MeshDensity',40)
title('primera derivada en x')
subplot(2,1,2)
fsurf(dfy,'MeshDensity',40)
title('primera derivada en y')
fyy=diff(f, y, 2); 
fxx=diff(f, x, 2);
fxy=diff(diff(f, x), y);
fyx=diff(diff(f, y), x);
a=5;
%% grafica curvas de nivel y gradiente
fx_n=matlabFunction(dfx,'Vars',{x , y});
fy_n=matlabFunction(dfy,'Vars',{x , y});
% Se define un mallado fino para las curvas de nivel
I = -5:.05:5; J = I; 
[X,Y] = meshgrid(I,J);
% Se define un mallado grueso para el gradiente
xx = -5:.2:5; yy = xx;
[XX,YY] = meshgrid(xx,yy);
% Se evalua la funcion y se calculan las 
% componentes del vector gradiente
Z = f_n(X,Y);
U = fx_n(XX,YY); 
V = fy_n(XX,YY);
% Se eligen los niveles (valores) 
% de las curvas a representar
levels = [-10:1:10]; 
figure
[c,h]=contour(X,Y,Z,levels);
clabel(c,h)
xlabel('x-axis')
ylabel('y-axis')
title('Curvas de nivel y campo gradiente')
hold on
quiver(XX,YY,U,V)
axis equal
%% Representacion del modulo del gradiente
modulo=sqrt(dfx^2+dfy^2);
x_min=-10;x_max=10;y_min=-2;y_max=2;
region_2D = [x_min x_max y_min y_max];
figure
fsurf(modulo,region_2D,'MeshDensity',40)
view(141,47)
title('Grafica del modulo del gradiente')
title('modulo del gradiente')
%% Grafica de la Hessiana por componentes
close all
figure
subplot(2,2,1)
fsurf(fxx,[-a,a,-a,a],'MeshDensity',40)
title('fxx: derivada parcial segunda en x')
xlabel('x');ylabel('y');
view(-41,46)
subplot(2,2,2)
fsurf(fxy,[-a,a,-a,a],'MeshDensity',40)
title('fxy: derivada parcial segunda cruzada')
xlabel('x');ylabel('y');
view(-71,64)
subplot(2,2,3)
fsurf(fyx,[-a,a,-a,a],'MeshDensity',40)
xlabel('x');ylabel('y');
title('fyx: derivada parcial segunda cruzada')
view(-71,64)
subplot(2,2,4)
fsurf(fyy,[-a,a,-a,a],'MeshDensity',40)
xlabel('x');ylabel('y');
title('fyy: derivada parcial segunda en y')
%%  campo de gradientes y de la hessiana
G1=[dfx,dfy];
H1=[fxx, fxy; fyx, fyy];
% definicion con los operadores gradient y hessian
G2=gradient(f,[x,y]);
H2=hessian(f,[x,y]);
% definicion con le operador jacobiano
G3= jacobian(f, [x, y]);
H3= jacobian(G3, [x, y]);
% Calculo norma y determinante 
norma=norm(G1);
determinante=fxx*fyy-fxy^2;
%% Representacion de la norma del gradiente 
% determinante de la hessiana
a=5;
figure
subplot(1,2,1)
fsurf(norma,[-a,a,-a,a],'MeshDensity',40)   %norma del determinante    
xlabel('x');ylabel('y');
title('norma del vector gradiente')
view(-41,46)
hold on
subplot(1,2,2)
fsurf(determinante,[-a,a,-a,a],'MeshDensity',40)    %determinante de la matriz hessiana
xlabel('x');ylabel('y');
title('determinante de la matriz hessiana')
view(46,45)
hold off
%% Restricciones de la funcion sobre la frontera 
f_S=f(x,-5);
f_N=f(x,5);
f_O=f(-5,y);
f_E=f(5,y);
figure
subplot(1,4,1)
ezplot(f_S,[-5,5])
title('frontera Sur')
xlabel('x')
subplot(1,4,2)
ezplot(f_N,[-5,5])
xlabel('x')
title('frontera Nord')
subplot(1,4,3)
ezplot(f_O,[-5,5])
ylabel('y')
title('frontera Oeste')
subplot(1,4,4)
ezplot(f_E,[-5,5])
ylabel('y')
title('frontera Este')
%% Aproximacion de Taylor  
%% Polinomio de Taylor de orden n centrado en [0,0]
a=5;
b=5;
n=2;    %el orden de la aproximacion es n-1  
T2=taylor(f,[x,y],'ExpansionPoint',[5 5],'Order',n);
I=[-a,a,-a,a];
J=[-b,b,-b,b];
fsurf(f,I,'MeshDensity', 40)
title('Superficie')
hold on
fsurf(T2,J,'FaceAlpha',0.5,'MeshDensity', 50)
hold off
title('Aproximacion lineal de Taylor')
%% Apartado Errores asbolutos y relativos en el punto de la aproximacion
valor_exacto=eval(subs(f,[x,y],P));
valor_lin=eval(subs(T2,[x,y],P));
err_abs_lin=abs(valor_exacto-valor_lin);
err_rel_lin=abs(valor_exacto-valor_lin)/abs(valor_exacto);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimizacion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Puntos criticos
%Hayo la derivadas 
fx=diff(f,x); fy=diff(f,y); 
[px,py]=solve(fx,fy); 
%% Clasificacion mediante el criterio
fxx(x,y)=diff(f,x,x);
H1=[fxx, fxy; fyx, fyy];
H(x,y)=hessian(f); 
minx=[];
maxx=[];
miny=[];
maxy=[];
psx=[];
psy=[];
for i=1:length(px)
if imag(px(i))==0 && imag(px(i))==0
H1=fxx(eval(px(i)),eval(py(i))); 
H2=det(H(eval(px(i)),eval(py(i))));
if(double(H1)>0 && double(H2)>0)
   minx=[minx double(px(i))];
   miny=[miny double(py(i))];
end
if(double(H1)<0 && double(H2)>0)
   maxx=[maxx double(px(i))];
   maxy=[maxy double(py(i))]; 
end
if(double(H2)<0)
   psx=[psx double(px(i))];
   psy=[psy double(px(i))];
end
end
end
%% Aproximacion lineal y cuadratica en los puntos criticos
disp('max')
fsurf(f,'MeshDensity',40) 
for i=1:length(maxx)
fprintf("(%d, %d)\n",maxx(i),maxy(i))
hold on
plot3(maxx(i),maxy(i),subs(f,[x,y],[maxx(i),maxy(i)]),'b*')
end
for i=1:length(minx)
fprintf("(%d, %d)\n",minx(i),miny(i))
hold on
plot3(minx(i),miny(i),subs(f,[x,y],[minx(i),miny(i)]),'r*')
end
for i=1:length(psx)
fprintf("(%d, %d)\n",psx(i),psy(i))
hold on
plot3(psx(i),psy(i),subs(f,[x,y],[psx(i),psy(i)]),'g*')
end
%% diferenciabilidad
%% Aproximacion lineal y cuadratica en los puntos criticos
%%
T2=zeros(6,1);
dir_x=y^2;
dir_y=x^2;
colx=zeros(10,1);
coly=zeros(10,1);
for i=1:length(px)
if imag(px(i))==0 
%{
aqui calculo diferenciabilidad de la funcion en los puntos maximos y
minimos utilizando el polinomio de taylor de orden 1 que equivale a el
plano tangente, de esta forma si su limite en el punto es 0 entonces
significara que es diferenciable, sin embargo estoy haciendo algun calculo
mal, por lo tanto dejo esta parte comentada
dife(x,y)=(f-(f(px(i),py(i))+fx*(x-px(i))+fy*(y-py(i))))/sqrt((x-px(i))^2+(y-py(i))^2);
L1=limit(limit(dife,x,px(i)),y,py(i))
L2=limit(limit(dife,y,px(i)),x,py(i))
if L2==0
%}   
% en caso afirmativo, meteremos el plano tangente en la lista    
T2(i)=taylor(f,[x,y],'ExpansionPoint',[px(i),py(i)],'Order',2);
colx(i)=px(i);
coly(i)=py(i);
end
end
%% superficie y de los planos tangentes
close all
fsurf(f)
for i=1:length(T2)
hold on
fsurf(T2(i),[colx(i)-1 colx(i)+1 coly(i)-1 coly(i)+1],'FaceAlpha',0.5,'MeshDensity', 3)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aparatado creativo
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
syms x y t
P_0=[maxx(1),maxy(1)];
P_f=[maxx(2),maxy(2)];
T_0=0;
T_f=5;
lambda=(t-T_0)/(T_f-T_0);
C=lambda*P_f +(1-lambda)*P_0;
long_2D=norm(P_f-P_0);
S=y*1.2*exp(-(x^2+y^2))*((-x^3+y^2)^3)*(y * x)^2;
x1=C(1); 
y1=C(2); 
z1=subs(S,[x,y], [x1,y1]);

r1=[maxx(1),maxy(1),z1];
X1=matlabFunction(x1);
Y1=matlabFunction(y1);
Z1=matlabFunction(z1);
figure
ezsurf(S)
shading interp
hold on
ezplot3(X1,Y1,Z1,[0,5,-20,20],'animate') % curva 3D
hold off