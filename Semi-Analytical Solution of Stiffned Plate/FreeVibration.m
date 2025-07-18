% Part 1 : Defination of Geometry
clc;clear;
%Case 1.1 single straight stiffener parallel to y axis
% a=0.5;ba=1;ha=0.01;b=a*ba;h=a*ha;
% bh=2;hh=4;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% Mx=8;Ny=8;
% load('.\BaseFunction\Case1.1\CCCC15');
% x1=0.4*a;y1=-b;x2=0.4*a;y2=0;x3=0.4*a;y3=b;ct=1;

%Case 1.2 single straight stiffener parallel to x axis
% a=0.5;ba=1;ha=0.02;b=a*ba;h=a*ha;
% bh=1;hh=1.5;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% Mx=8;Ny=8;
% load('.\BaseFunction\Case1.2\SSSS15');
% x1=-a;y1=0;x2=0;y2=0;x3=a;y3=0;ct=1;

%Case 1.3 cross stiffeners
% a=0.5;ba=1;ha=0.02;b=a*ba;h=a*ha;
% bh=1;hh=1.5;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% Mx=8;Ny=8;
% load('.\BaseFunction\Case1.2\SSSS15');
% x1=[-a,0];y1=[0,-b];x2=[0,0];y2=[0,0];x3=[a,0];y3=[0,b];ct=[1,1];

%Case 1.4 two stiffeners along x aixs and two stiffeners alon y axis
% a=0.5;ba=0.5;ha=0.01;b=a*ba;h=a*ha;
% bh=1;hh=1;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% Mx=8;Ny=8;
% load('.\BaseFunction\Case1.4\SSSS10');
% x1=[-a,-a,-a/3,a/3];y1=[-b/3,b/3,-b,-b];
% x2=[0,0,-a/3,a/3];y2=[-b/3,b/3,0,0];
% x3=[a,a,-a/3,a/3];y3=[-b/3,b/3,b,b];ct=[1,1,1,1];

%Case 2.1 skew cross stiffeners
% a=0.5;ba=1;ha=0.01;b=a*ba;h=a*ha;
% bh=2;hh=4;bs=h*bh;hs=h*hh;e=0;Js=(bs*hs^3+bs^3*hs)/12;
% E=70E9;v=0.3;rho=2700;
% MN=15;Mx=15;Ny=15;
% load('.\BaseFunction\Case2.1\SSSS15');
% x1=[-0.5,-0.5];y1=[-0.5,0.5];
% x2=[0,0];y2=[0,0];
% x3=[0.5,0.5];y3=[0.5,-0.5];ct=[1,1];

% Case 2.2 tow skew stiffeners
% a=0.2032/2;ba=1;ha=0.001/a;b=a*ba;h=a*ha;
% bh=2.032;hh=1.5;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case2.1\SSSS15');
% x1=[0-a,0.05-a];y1=[0.1-b,0.2032-b];
% x2=[0.1016-a,0.1-a];y2=[0.125-b,0.1016-b];
% x3=[0.2032-a,0.15-a];y3=[0.15-b,0-b];ct=[1,1];

% Case 3.1 single parabolic stiffener 
% a=0.6069/2;ba=0.7112/0.6069;ha=0.00625/a;b=a*ba;h=a*ha;
% bh=4.06/6.25;hh=22.4/6.25;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case3.1\SSSS15');
% x1=-a;y1=1.2*(0.6-1.5*0)^2+0.05-b;
% x2=0;y2=1.2*(0.6-1.5*a)^2+0.05-b;
% x3=a;y3=1.2*(0.6-1.5*a*2)^2+0.05-b;ct=1;

% Case 3.2 single parabolic stiffener 
% a=0.6069/2;ba=0.7112/0.6069;ha=0.00625/a;b=a*ba;h=a*ha;
% bh=4.06/6.25;hh=22.4/6.25;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case3.2\SSSS15');
% x1=-a;y1=-30*(0.1-0.35*0)^2+0.4-b;
% x2=0;y2=-30*(0.1-0.35*a)^2+0.4-b;
% x3=a;y3=-30*(0.1-0.35*a*2)^2+0.4-b;ct=1;

% Case 3.3 single parabolic stiffener 
% a=0.6069/2;ba=0.7112/0.6069;ha=0.00625/a;b=a*ba;h=a*ha;
% bh=4.06/6.25;hh=4;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case3.3\SSSS15');
% x1=-a;y1=0.4820-b;x2=0.3034-a;y2=0.0752-b;x3=0.6069-a;y3=0.1656-b;ct=2;

%Case 3.4 single parabolic stiffener 
% a=0.7;ba=0.5;ha=0.00625/a;b=a*ba;h=a*ha;
% bh=4.06/6.25;hh=4;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case3.4\CCCC20');
% x1=0-a;y1=0.482-b;x2=0.7-a;y2=0.6752-b;x3=1.4-a;y3=0.1656-b;ct=2;


% Case 3.5 two parabolic stiffener 
a=0.6;ba=7/12;ha=0.00625/a;b=a*ba;h=a*ha;
bh=4.06/6.25;hh=3;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
E=73E9;v=0.33;rho=2837;
Mx=15;Ny=15;
load('.\BaseFunction\Case3.5\CCCC20');
x1=[0-a,0-a];y1=[0.482-b,0.234-b];ct=2;
x2=[0.6-a,0.6-a];y2=[0.6752-b,0.081-b];
x3=[1.2-a,1.2-a];y3=[0.1356-b,0.623-b];ct=[2,2];

%Case 3.6 two parabolic stiffener 
% a=0.6;ba=2/3;ha=0.00625/a;b=a*ba;h=a*ha;
% bh=4.06/6.25;hh=4;bs=h*bh;hs=h*hh;e=(h+hs)/2;Js=bs^3*hs/3;
% E=73E9;v=0.33;rho=2837;
% Mx=15;Ny=15;
% load('.\BaseFunction\Case3.5\CCCC20');
% x1=[0-a,0.4-a];y1=[0.482-b,0-b];
% x2=[0.6-a,0.8-a];y2=[0.6752-b,0.2666-b];
% x3=[1.2-a,0.8-a];y3=[0.1356-b,0.8-b];ct=[2,2];
%%
% Part 2 ： Sovling 
Dimp=[a,ba,ha];Dims=[bh,hh,e,Js];Mat=[E,v,rho];
Ip=[x1;y1;x2;y2;x3;y3;ct];
load('.\BaseFunction\GL50');% Gauss-Lobatto integral points and coefficients
PS=parameter(Dimp,Dims,Mat);
[omegas0,Xeigs,As,Yeigs,Bs]=selectmodes(omegas,Xeigs,As,Yeigs,Bs,Mx,Ny);
[Mp,Kp]=Pmatrix(PS,Xeigs,As,Yeigs,Bs,xi,C);
[Ms,Ks]=Smatrix(PS,Xeigs,As,Yeigs,Bs,Ip,@parabola,xi,C);
M=Mp+Ms;K=Kp+Ks;
[Modes,Omegas2]=eig(K,M,'qz');
[omegas,v]=sort(sqrt(diag(Omegas2)));
modes=Modes(:,v);
%%
% Part 3 ： Showing frequency and mode
D=PS(9);k1=a^2*sqrt(rho*h/D);

%kk=1:12;% kk defines the desired mode
kk=[1:6,10,20,30,50];

%omegas(kk)/k1/pi/2 % Hz
omegas(kk)*4*ba^2/pi^2 % dimensionless

qs=kk;Figy=3;Figx=4;Fign=1;Ft=2;
%Ft==1 represents w and Ft==other value represents abs(w)
Ws=plotmodes(PS,Xeigs,As,Yeigs,Bs,modes,qs,Fign,Figx,Figy,Ft);