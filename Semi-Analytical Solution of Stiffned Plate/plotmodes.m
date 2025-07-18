function Ws=plotmodes(PS,Xeigs,As,Yeigs,Bs,modes,qs,Fign,Figx,Figy,Ft)
%plot the modes for stiffened panel
%input : 
% PS  the parameters of the panel, stiffeners and airflow
% Xeigs,Yeigs  eigenvalue of phi(x) and varphi(y)
% As,Bs  the  mode coefficients phi(x) and varphi(y)
% modes  the amplitude of generalized coordinates
% qs defines the modes that need to be drawn
% Fign defines figure
% Figx,Figy defines the dimensions of subplots
x=(-1:0.02:1)';y=x;a=PS(1);b=PS(2)*a;a2=a*x;b2=b*y;
m=length(qs);[~,n]=size(Xeigs);k=length(x);
Ws=cell(m,1);
f1=figure(Fign);
for i=1:m
    mode=modes(:,qs(i));
    W=zeros(k,k);
    for j=1:n
        phix=As(1,j).*cosh(Xeigs(1,j).*x)+As(2,j).*sinh(Xeigs(1,j).*x)+...
            As(3,j).*cosh(Xeigs(2,j).*x)+As(4,j).*sinh(Xeigs(2,j).*x);
        phiy=Bs(1,j).*cosh(Yeigs(1,j).*y)+Bs(2,j).*sinh(Yeigs(1,j).*y)+...
            Bs(3,j).*cosh(Yeigs(2,j).*y)+Bs(4,j).*sinh(Yeigs(2,j).*y);
        W=mode(j)*phiy*phix'+W;
    end
    Ws{i,1}=W./max(max(abs(W)));
    subplot(Figy,Figx,i)
    [X,Y]=meshgrid(a2,b2);
    if Ft==1
        surf(X,Y,0.02*real(W));
    else
        surf(X,Y,0.02*abs(real(W)));
    end
    axis vis3d;axis equal;grid on;axis off;
    colormap jet;shading interp;f1.Color = 'white';
    view(2);
end
end