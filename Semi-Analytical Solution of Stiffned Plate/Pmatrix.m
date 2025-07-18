function [MP,KP]=Pmatrix(PS,Xeigs,As,Yeigs,Bs,xi,C)
% get the mass matrix Mp and stiffness matrix Kp for the plate
% input :
% PS the parameters of the panel, stiffeners and airflow
% Xeigs,Yeigs  eigenvalues of trial functions phi(x) and varphi(y)
% As,Bs  the  mode coefficients of trial functions phi(x) and varphi(y)
% xi the Gauss-Lobatto integral points (-1,1)
% C  the Gauss-Lobatto integral coefficients

M=length(xi);[~,N]=size(Xeigs);
% M is the number of Gauss-Lobatto integral points
% N is the number of base functions
phix=zeros(M,N);phiy=zeros(M,N);dphix=zeros(M,N);dphiy=zeros(M,N);
ddphix=zeros(M,N);ddphiy=zeros(M,N);
for i=1:N
    A=As(:,i);B=Bs(:,i);Xeig=Xeigs(:,i);Yeig=Yeigs(:,i);
    [phixi,phiyi,dphixi,dphiyi,ddphixi,ddphiyi]=dphixy(A,Xeig,B,Yeig,xi);
    phix(:,i)=phixi;phiy(:,i)=phiyi;
    dphix(:,i)=dphixi;dphiy(:,i)=dphiyi;
    ddphix(:,i)=ddphixi;ddphiy(:,i)=ddphiyi;
end
%
ba=PS(2);v=PS(12);
k1=1/ba^4;k2=v/ba^2;k3=2*(1-v)/ba^2;
KP=zeros(N,N);MP=zeros(N,N);
for i=1:N
    for j=i:N
        px00=sum(phix(:,i).*phix(:,j).*C);
        py00=sum(phiy(:,i).*phiy(:,j).*C);
        
        px11=sum(dphix(:,i).*dphix(:,j).*C);
        py11=sum(dphiy(:,i).*dphiy(:,j).*C);
        
        px22=sum(ddphix(:,i).*ddphix(:,j).*C);
        py22=sum(ddphiy(:,i).*ddphiy(:,j).*C);
        
        px20=sum(ddphix(:,i).*phix(:,j).*C);
        py20=sum(ddphiy(:,i).*phiy(:,j).*C);
        
        px02=sum(phix(:,i).*ddphix(:,j).*C);
        py02=sum(phiy(:,i).*ddphiy(:,j).*C);
        
        kp=px22*py00+k1*px00*py22+k2*(px20*py02+px02*py20)+k3*px11*py11;
        mp=px00*py00;
        KP(i,j)=kp;MP(i,j)=mp;
        MP(j,i)=MP(i,j);KP(j,i)=KP(i,j);
    end
end
end

function [phix,phiy,dphix,dphiy,ddphix,ddphiy]=dphixy(A,Xeig,B,Yeig,xi)
% get the derivative of trial function w=phix*phiy
x=xi;y=xi;

phix=A(1)*cosh(Xeig(1).*x)+A(2)*sinh(Xeig(1).*x)+...
     A(3)*cosh(Xeig(2).*x)+A(4)*sinh(Xeig(2).*x);
phiy=B(1).*cosh(Yeig(1).*y)+B(2).*sinh(Yeig(1).*y)+...
     B(3).*cosh(Yeig(2).*y)+B(4).*sinh(Yeig(2).*y);
phix=real(phix);phiy=real(phiy);
 
dphix=A(1)*Xeig(1)*sinh(Xeig(1).*x)+A(2)*Xeig(1)*cosh(Xeig(1).*x)+...
      A(3)*Xeig(2)*sinh(Xeig(2).*x)+A(4)*Xeig(2)*cosh(Xeig(2).*x);
dphiy=B(1)*Yeig(1)*sinh(Yeig(1).*y)+B(2)*Yeig(1)*cosh(Yeig(1).*y)+...
      B(3)*Yeig(2)*sinh(Yeig(2).*y)+B(4)*Yeig(2)*cosh(Yeig(2).*y);  
dphix=real(dphix);dphiy=real(dphiy);
% A,Xeig,B,Yeig are obtained in the dimensionless coordinates x and y

ddphix=A(1)*(Xeig(1))^2*cosh(Xeig(1).*x)+A(2)*(Xeig(1))^2*sinh(Xeig(1).*x)+...
       A(3)*(Xeig(2))^2*cosh(Xeig(2).*x)+A(4)*(Xeig(2))^2*sinh(Xeig(2).*x);
ddphiy=B(1)*(Yeig(1))^2*cosh(Yeig(1).*y)+B(2)*(Yeig(1))^2*sinh(Yeig(1).*y)+...
       B(3)*(Yeig(2))^2*cosh(Yeig(2).*y)+B(4)*(Yeig(2))^2*sinh(Yeig(2).*y);
ddphix=real(ddphix);ddphiy=real(ddphiy);
end


