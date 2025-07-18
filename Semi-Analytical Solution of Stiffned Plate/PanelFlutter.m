% Part 1 : Defination of Geometry
clear;clc;

%Case f.1 single straight stiffener parallel to y axis
% a=0.5;ba=1;ha=0.01;b=a*ba;h=a*ha;
% bh=2;hh=4;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% Mx=8;Ny=8;
% load('.\BaseFunction\Casef.1\SSSS15');
% x1=0.4*a;y1=-b;x2=0.4*a;y2=0;x3=0.4*a;y3=b;

%Case f.2  straight stiffeners
a=0.5;ba=1;ha=0.05;bh=0.4;b=a*ba;h=a*ha;
hh=1.5;b=a*ba;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
E=70E9;v=0.3;rho=2700;
Mx=10;Ny=10;
load('.\BaseFunction\Casef.2\CCCC15');
%x1=-a;y1=0;x2=0;y2=0;x3=a;y3=0;ct=1;% one central x stiffener
%x1=0;y1=-b;x2=0;y2=0;x3=0;y3=b;ct=1;% one central y stiffener
%x1=[-a,-a];y1=[-b/3,b/3];x2=[0,0];y2=[-b/3,b/3];x3=[a,a];y3=[-b/3,b/3];ct=[1,1];% two x stiffeners
%x1=[-a/3,a/3];y1=[-b,-b];x2=[-a/3,a/3];y2=[0,0];x3=[-a/3,a/3];y3=[b,b];ct=[1,1];% two y stiffeners
x1=[-a,0];y1=[0,-b];x2=[0,0];y2=[0,0];x3=[a,0];y3=[0,b];ct=[1,1];% cross stiffeners
%%
% Part 2 ： Sovling 
Dimp=[a,ba,ha];Dims=[bh,hh,e,Js];Mat=[E,v,rho];
Ip=[x1;y1;x2;y2;x3;y3;ct];
load('.\BaseFunction\GL50');% Gauss-Lobatto integral points and coefficients
PS=parameter(Dimp,Dims,Mat);
[omegas0,Xeigs,As,Yeigs,Bs]=selectmodes(omegas,Xeigs,As,Yeigs,Bs,Mx,Ny);
[Mp,Kp]=Pmatrix(PS,Xeigs,As,Yeigs,Bs,xi,C);
[Ms,Ks]=Smatrix(PS,Xeigs,As,Yeigs,Bs,Ip,@parabola,xi,C);
Ka=Amatrix(Xeigs,As,Yeigs,Bs,xi,C);
gamma0=0:0.1:200;
Omegas0=EigSolve(Mp,Kp,Ms,Ks,Ka,gamma0);
%%
% Part 3 ： Showing flutter frequency
Num=(Mx*Ny+1):(Mx*Ny+4);
%plotfrequency(gamma0,Omegas0/k1/pi/2,Nplot1,3);
plotfrequency(gamma0*8,Omegas0*4,Num,3);