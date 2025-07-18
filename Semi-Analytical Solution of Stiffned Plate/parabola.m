function [x,y,dx,dy]=parabola(curves,xi)
% get the vlaue, the first derivative and the second derivative of
% a parabola (or a skew line) x=x(xi),y=y(xi) for given xi
% input : 
% curves=[x1;y1;x2;y2;x3;y3;ct];
% (x1,y2),(x2,y2),(x3,y3) define the parabolic curve, ct==1 denotes 
% interpolated points and ct==2 denotes control points.
% xi the Gaussian-Lobatto integral point
m=length(xi);[~,n]=size(curves);
x=zeros(m,n);y=zeros(m,n);
dx=zeros(m,n);dy=zeros(m,n);
for k=1:n
if curves(end,k)==1
    N1=0.5*(xi-1).*xi;N2=(1-xi).*(1+xi);N3=0.5.*(1+xi).*xi;
    dN1=xi-0.5;dN2=-2.*xi;dN3=xi+0.5;
else
    N1=0.25*(xi-1).^2;N2=0.5*(1-xi).*(1+xi);N3=0.25.*(1+xi).^2;
    dN1=0.5*(xi-1);dN2=-xi;dN3=0.5*(xi+1);
end
    x1=curves(1,k);y1=curves(2,k);
    x2=curves(3,k);y2=curves(4,k);
    x3=curves(5,k);y3=curves(6,k);
    x(:,k)=N1.*x1+N2.*x2+N3.*x3;
    y(:,k)=N1.*y1+N2.*y2+N3.*y3;
    dx(:,k)=dN1.*x1+dN2.*x2+dN3.*x3;
    dy(:,k)=dN1.*y1+dN2.*y2+dN3.*y3;
end
end