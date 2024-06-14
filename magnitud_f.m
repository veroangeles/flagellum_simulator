function  magnitud_f
%clc
a=0.1;
e=(3/2)*a;
x01=0;
x02=0;
x03=0;
v=[0; 0; 1];
s=stoke(e,x01,x02,x03,x01,x02,x03)
f=s\v
err=abs(f(3)-6*pi*a)
end