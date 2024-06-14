function coupling 
% a es el radio de la esfera
a=0.1;
%e es el parametro de regularizacion para que la 
%el esokelet regularizado de valores cercanos a los de la esfera
e=(3/2)*a;
% x1, x2 son las posiciones de los centros de las dos esferas.
x1=[0;0;0];
x2=[0;1;0];
% velocidades de las dos esferas.
U1=[0;0;1];
U2=[0;0;1];
% la funcion Gij=stoke(e,xi(1),xi(2),xi(3),xj(1),xj(2),xj(3)) calcula el stokelet
% quenerado por una singularidad en xj  en xi
G11=stoke(e,x1(1),x1(2),x1(3),x1(1),x1(2),x1(3));
G12=stoke(e,x1(1),x1(2),x1(3),x2(1),x2(2),x2(3));
G21=stoke(e,x2(1),x2(2),x2(3),x1(1),x1(2),x1(3));
G22=stoke(e,x2(1),x2(2),x2(3),x2(1),x2(2),x2(3));
%creamos la matriz G para calcular las magnitudes de las dos singularidades
%en los puntos x1 y x2. Las primeras tres entradas de F son la magnitud de
%la primera sigularidad y las otras tres entradas son la magnitud de la
%segunda
G=[G11 G12; G21 G22];
U = [U1;U2];

%F= inv(G)*U
F = G\U
%size(U)
%size(F)
end



