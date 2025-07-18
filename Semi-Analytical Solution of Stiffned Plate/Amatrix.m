function KA=Amatrix(Xeigs,As,Yeigs,Bs,xi,C)
% get the areodynamic stiffness matrix KA
% input: 
% PS the parameters of the panel, stiffeners and airflow
% Xeig,Yeig  eigenvalue of phi(x) and varphi(y)
% A,B the  mode coefficients phi(x) and varphi(y)
% xi  the Gauss-Lobatto integral points (-1,1)
% C   the Gauss-Lobatto integral coefficients
M=length(xi);[~,N]=size(Xeigs);
% M is the number of Gauss-Lobatto integral points
% N is the number of base functions
phix=zeros(M,N);phiy=zeros(M,N);dphix=zeros(M,N);dphiy=zeros(M,N);
for i=1:N
    A=As(:,i);B=Bs(:,i);Xeig=Xeigs(:,i);Yeig=Yeigs(:,i);
    [phixi,phiyi,dphixi,dphiyi]=dphixy(A,Xeig,B,Yeig,xi);
    phix(:,i)=phixi;phiy(:,i)=phiyi;
    dphix(:,i)=dphixi;dphiy(:,i)=dphiyi;
end

KA=zeros(N,N);
for i=1:N
    for j=1:N
        py00=sum(phiy(:,i).*phiy(:,j).*C);      
        px10=sum(dphix(:,j).*phix(:,i).*C);
        KA(i,j)=px10*py00;
    end
end
end

function [phix,phiy,dphix,dphiy]=dphixy(A,Xeig,B,Yeig,xi)
% get the derivative of base function w=phix*phiy
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
% A,Xeig,B,Yeig are obtained by the dimensionless coordinates x and y
end


