function  plotfrequency(gamma,Omegas0,Num,Fign)
%display the relation between Omega and gamma
%input : gamma the arodynamic stiffness coefficient
%        Omegas0 the flutter eigenvalue
%        Num determines the number of Omega that is required to display
%        Fign define the figure
m=length(Num);
figure(Fign)
hold on
for i=1:m
    plot(gamma,real(Omegas0(:,Num(i))))
end
hold off
grid on;

figure(Fign+1)
hold on
for i=1:m
    plot(gamma,imag(Omegas0(:,Num(i))))
end
hold off
grid on;

    