function [MS,KS]=Smatrix(PS,Xeigs,As,Yeigs,Bs,Ip,curve,xi,C)
% get the stiffness and mass matrix for stiffeners
% input:
% PS the parameters of the panel, stiffeners and airflow
% Xeigs,Yeigs,As,Bs define the base functions
% Ip the interpolation or control points of stiffener
% curve a function the defines the geometry of stiffeners
% xi the Gauss-Lobatto integral point
% C  the Gauss-Lobatto integral coefficients

[~,S]=size(Ip);[~,N]=size(Xeigs);
% S is the number of stiffener, N is the number of trial functions
[Xs,Ys,dXs,dYs]=curve(Ip,xi);
a=PS(1);ba=PS(2);ha=PS(3);
A=PS(4);Ae2=PS(5);In=PS(6);Ib=PS(7);Js=PS(8);
D=PS(9);Es=PS(10);Gs=PS(11);
k0=a^2/(ba*D);m0=1/(a.^3*ba*ha);
kn=Es*(In+Ae2)*k0;kl=Gs*Js*k0;
mw=A*m0;ml=(In+Ib)*m0;
MS=zeros(N,N);KS=zeros(N,N);
for i=1:S
    X=Xs(:,i);Y=Ys(:,i);dX=dXs(:,i);dY=dYs(:,i);
    [Kn,Kl,Mw,Ml]=smatrix(PS,Xeigs,As,Yeigs,Bs,X,Y,dX,dY,C);
    KS=KS+kn*Kn+kl*Kl;
    MS=MS+mw*Mw+ml*Ml;
end
end

function [Kn,Kl,Mw,Ml]=smatrix(PS,Xeigs,As,Yeigs,Bs,x,y,dx,dy,C)
% get the stiffness and mass matrix for single stiffener
[~,N]=size(Xeigs);M=length(x);
J=sqrt(dx.^2+dy.^2);lx=dx./J;ly=dy./J;
% N is the number of base functions, M is the number of Gauss-Lobatto integral point
w=zeros(M,N);wl=zeros(M,N);wn=zeros(M,N);wlx=zeros(M,N);wnx=zeros(M,N);
for i=1:N
    A=As(:,i);B=Bs(:,i);Xeig=Xeigs(:,i);Yeig=Yeigs(:,i);
    [wi,wxi,wyi,wxxi,wyyi,wxyi]=dW(PS,A,Xeig,B,Yeig,x,y);
    w(:,i)=wi;
    wl(:,i)=lx.*wxi+ly.*wyi;
    wn(:,i)=lx.*wyi-ly.*wxi;
    wlx(:,i)=lx.^2.*wxxi+ly.^2.*wyyi+2.*lx.*ly.*wxyi;
    wnx(:,i)=lx.*ly.*(wyyi-wxxi)+(lx.^2-ly.^2).*wxyi;
end
%
Kn=zeros(N,N);Kl=zeros(N,N);Mw=zeros(N,N);Ml=zeros(N,N);
for i=1:N
    for j=i:N
        chinij=wlx(:,i).*wlx(:,j);chilij=wnx(:,i).*wnx(:,j);
        wij=w(:,i).*w(:,j);wnij=wn(:,i).*wn(:,j);
        Kn(i,j)=sum(chinij.*J.*C);Kl(i,j)=sum(chilij.*J.*C);
        Mw(i,j)=sum(wij.*J.*C);Ml(i,j)=sum(wnij.*J.*C);
        if j>i
            Kn(j,i)=Kn(i,j);Kl(j,i)=Kl(i,j);
            Mw(j,i)=Mw(i,j);Ml(j,i)=Ml(i,j);
        end
    end
end
end


function [w,wx,wy,wxx,wyy,wxy]=dW(PS,A,Xeig,B,Yeig,x0,y0)
% get the derivative of trial function w=phix*phiy
a=PS(1);b=a*PS(2);x=x0./a;y=y0./b;
% here, the symbol 'x' and 'y' represent the dimensionless coordinate 
% x=x0/a and y=y0/b.the symbol 'x0' and 'y0' represent the
% phisical (dimensional) coordinate
phix=A(1)*cosh(Xeig(1).*x)+A(2)*sinh(Xeig(1).*x)+...
     A(3)*cosh(Xeig(2).*x)+A(4)*sinh(Xeig(2).*x);
phiy=B(1).*cosh(Yeig(1).*y)+B(2).*sinh(Yeig(1).*y)+...
     B(3).*cosh(Yeig(2).*y)+B(4).*sinh(Yeig(2).*y);
phix=real(phix);phiy=real(phiy);
 
dphix=A(1)*Xeig(1)*sinh(Xeig(1).*x)+A(2)*Xeig(1)*cosh(Xeig(1).*x)+...
      A(3)*Xeig(2)*sinh(Xeig(2).*x)+A(4)*Xeig(2)*cosh(Xeig(2).*x);
dphiy=B(1)*Yeig(1)*sinh(Yeig(1).*y)+B(2)*Yeig(1)*cosh(Yeig(1).*y)+...
      B(3)*Yeig(2)*sinh(Yeig(2).*y)+B(4)*Yeig(2)*cosh(Yeig(2).*y);  
dphix=real(dphix)./a;dphiy=real(dphiy)./b;
% A,Xeig,B,Yeig are obtained in the dimensionless coordinates x and y

ddphix=A(1)*(Xeig(1))^2*cosh(Xeig(1).*x)+A(2)*(Xeig(1))^2*sinh(Xeig(1).*x)+...
       A(3)*(Xeig(2))^2*cosh(Xeig(2).*x)+A(4)*(Xeig(2))^2*sinh(Xeig(2).*x);
ddphiy=B(1)*(Yeig(1))^2*cosh(Yeig(1).*y)+B(2)*(Yeig(1))^2*sinh(Yeig(1).*y)+...
       B(3)*(Yeig(2))^2*cosh(Yeig(2).*y)+B(4)*(Yeig(2))^2*sinh(Yeig(2).*y);
ddphix=real(ddphix)./a^2;ddphiy=real(ddphiy)./b^2;

w=phix.*phiy;wx=dphix.*phiy;wy=phix.*dphiy;
wxx=ddphix.*phiy;wyy=phix.*ddphiy;wxy=dphix.*dphiy;
end



