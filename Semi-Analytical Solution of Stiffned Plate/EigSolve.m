function Omegas=EigSolve(Mp,Kp,Ms,Ks,Ka,gamma0)
%slove the eigenvalue equation for stiffened panel flutter
%input : 
% PS the parameters of panel and material
% Mp,Ms   the mass matrixes of plate and stiffener
% Kp,Ks   the stiffness matrix of plate and stiffener
% Ka      the aerodynamic the stiffness matrix
% gamma0   the dimensionless coefficient of aerodynamic stiffness
%output : 
% Omegas the flutter eigenvalue
% Modes  the flutter modes
% Num    recording the location of modes after sorting the frquency
n=length(gamma0);[~,m]=size(Mp);
M=Mp+Ms;Kps=Kp+Ks;
gamma=gamma0;alpha=0.1*sqrt(gamma0/2);
Omegas=zeros(n,2*m);
for i=1:n
    K=Kps+gamma(i)*Ka;R=alpha(i).*Mp;
    gamma(i)
    I1=eye(m);I0=zeros(m);
    Msum=[M,I0;I0,I1];Ksum=[R,K;-I1,I0];
    [~,Omega]=eig(-Ksum,Msum,'qz');
    Omegas(i,:)=conj((diag(Omega))');
end
[Gomega,~]=sort(-Omegas.*1i,2,'ComparisonMethod','real');
Omegas=1i.*Gomega;
end