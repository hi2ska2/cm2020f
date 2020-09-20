q = 1.602192e-19; 
eps0 = 8.854187e-12; 
Delta_x = 0.1e-9;
N = 67; 
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e23; % 1e18/cm^3, we can change it to 1e23,24,25 
interface1 = 9;
interface2 = 59;

A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
for i=2:N-1
     if     (i < interface1) 
         A(i,i-1) = eps_ox;    A(i,i) = -2*eps_ox;         A(i,i+1) = eps_ox;
     end
     if (i == interface1) 
         A(i,i-1) = eps_ox;    A(i,i) = -eps_ox-eps_si;    A(i,i+1) = eps_si;
     end
     if (i < interface2)   
         A(i,i-1) = eps_si;    A(i,i) = -2*eps_si;         A(i,i+1) = eps_si;
     end
     if (i == interface2) 
         A(i,i-1) = eps_si;    A(i,i) = -eps_si-eps_ox;    A(i,i+1) = eps_ox;
     end
     if (i > interface2) 
         A(i,i-1) = eps_ox;    A(i,i) = -2*eps_ox;         A(i,i+1) = eps_ox;
     end
end
b= zeros(N,1);
for i=interface1:interface2
    if (i==interface1) b(i,1)= Delta_x*Delta_x*q*Nacc/eps0*0.5;
    elseif (i==interface2) b(i,1)= Delta_x*Delta_x*q*Nacc/eps0*0.5;
    else b(i,1) = Delta_x*Delta_x*q*Nacc/eps0;
        end
end

x=Delta_x*[0:(N-1)];
y=inv(A)*b;
plot(x,y);
xlabel('x[nm]')
ylabel('V')

y_analytic = zeros(N,1);
%analytic solution
for i=1:interface1
    y_analytic(i,1)=1/2*x(i)*-3*0.8e-9*q*Nacc*5e-9/(eps0*eps_si)/0.8e-9;
end
for i= interface1+1 : interface2-1
    y_analytic(i,1)=1/2*Nacc*(x(i)-0.8e-9)*(x(i)-5.8e-9)/(eps0*eps_si)+1/2*-3*0.8e-9*q*Nacc*5e-9/(eps0*eps_si);
end
for i= interface2 : N
    y_analytic(i,1)=1/2*(x(i)-1.6e-9-5e-9)*-3*0.8e-9*q*Nacc*5e-9/(eps0*eps_si)/0.8e-9;
end
