function  [omega,Xeig,A,Yeig,B]=selectmodes(omegas,Xeigs,As,Yeigs,Bs,Mx,Ny)
%selecting modes as the base function
% omegas,omega the natural frequency
% Xeigs,Yeigs,Xeig,Yeig the eigenvalue of natural mode
% As,Bs,A,B the coefficients of natural mode
% Mx,Ny the number of nodal line along x,y direction, respectively
k=Mx*Ny;omega=zeros(k,1);Xeig=zeros(2,k);Yeig=zeros(2,k);
A=zeros(4,k);B=zeros(4,k);
[~,K]=size(Xeigs);Ns=sqrt(K);
for i=1:Ny
    for j=1:Mx
        omega((i-1)*Mx+j)=omegas((i-1)*Ns+j,1);
        Xeig(:,(i-1)*Mx+j)=Xeigs(:,(i-1)*Ns+j);
        Yeig(:,(i-1)*Mx+j)=Yeigs(:,(i-1)*Ns+j);
        A(:,(i-1)*Mx+j)=As(:,(i-1)*Ns+j);
        B(:,(i-1)*Mx+j)=Bs(:,(i-1)*Ns+j);
    end
end