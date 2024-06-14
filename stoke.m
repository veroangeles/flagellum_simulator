function[g]=stoke(e,x1,x2,x3,x01,x02,x03)
r=sqrt((x1-x01)^2+(x2-x02)^2+(x3-x03)^2);
%e=0.1;
re=sqrt(r^2+e^2);
%Dado que la matriz es simetrica solo obtengo 6 entradas 
mu=100;
g11=(r^2+2*e^2+(x1-x01)*(x1-x01))/(mu*8*pi*re^(3));
g12=((x1-x01)*(x2-x02))/(mu*8*pi*re^(3));
g13=((x1-x01)*(x3-x03))/(mu*8*pi*re^(3));
%g21=((x2-x02)*(x1-x01))/(8*pi*re^3);
g22=(r^2+2*e^2+(x2-x02)*(x2-x02))/(mu*8*pi*re^(3));
g23=((x2-x02)*(x3-x03))/(8*pi*re^(3));
%g31=((x3-x03)*(x1-x01))/(8*pi*re^3);
%g32=((x3-x03)*(x2-x02))/(8*pi*re^3);
g33=(r^2+2*e^2+(x3-x03)*(x3-x03))/(mu*8*pi*re^(3));
g=[g11 g12 g13; g12 g22 g23; g13 g23 g33];
end