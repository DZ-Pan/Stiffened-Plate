function PS=parameter(Dimp,Dims,Mat)
%define the parameters of the plate, stiffener and airflow
%input  : 
% Dimp=[a,ba,ha]      the parameters of plate
% Dims=[bh,hh,e,Js]   the parameters of stiffener
% Mat=[E,v,rho]       the paramters of material
% a  the half length of the plate
% ba the ratio of plate width and plate length
% ha the ratio of plate height and plate length
% bh the ration of beam width and plate height
% hh the ration of beam height and plate height
% e  the eccentricity of stiffener
% Js the rotational stiffness of stiffener

a=Dimp(1);ba=Dimp(2);ha=Dimp(3);
bh=Dims(1);hh=Dims(2);e=Dims(3);Js=Dims(4);
h=a*ha;bs=h*bh;hs=h*hh;
E=Mat(1);v=Mat(2);rho=Mat(3);G=E/(2*(1+v));
D=E*h^3/(12*(1-v^2));
As=bs*hs;Ae2=bs*hs*e^2;In=bs*hs^3/12;Ib=bs^3*hs/12;
%parameters of airflow
c=340;rhoa=1.205;alpha=a^2*rhoa*c/sqrt(rho*h*D);

PS=[a;ba;ha;As;Ae2;In;Ib;Js;D;E;G;v;alpha];
end