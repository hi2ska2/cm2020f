clear;
q = 1.602192e-19; 
eps0 = 8.854187e-12; 
Delta_x = 0.1e-9;
N = 67; 
T=300.0;
k = 1.380662e-23;
ni = 1.075e16;
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e24; % 1e18/cm^3, we can change it to 1e23,24,25 
interface1 = 9;
interface2 = 59;
Applied_Voltage = 0:0.1:1;
A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;

for i=2:(N-1)
     if     (i < interface1) 
         A(i,i-1) = eps_ox;    A(i,i) = -2*eps_ox;         A(i,i+1) = eps_ox;
     elseif (i == interface1) 
         A(i,i-1) = eps_ox;    A(i,i) = -eps_ox-eps_si;    A(i,i+1) = eps_si;
     elseif (i < interface2)   
         A(i,i-1) = eps_si;    A(i,i) = -2*eps_si;         A(i,i+1) = eps_si;
     elseif (i == interface2) 
         A(i,i-1) = eps_si;    A(i,i) = -eps_si-eps_ox;    A(i,i+1) = eps_ox;
     elseif (i > interface2) 
         A(i,i-1) = eps_ox;    A(i,i) = -2*eps_ox;         A(i,i+1) = eps_ox;
     end
end
b= zeros(N,11);
for j=1:11
b(1,j)=0.33374-(j-1)/10; 
b(N,j)=0.33374-(j-1)/10;
for i=interface1:interface2
    if (i==interface1) b(i,j)= Delta_x*Delta_x*q*Nacc/eps0*0.5;
    elseif (i==interface2) b(i,j)=  Delta_x*Delta_x*q*Nacc/eps0*0.5;
    else b(i,j) = Delta_x*Delta_x*q*Nacc/eps0;
        end
end
phi(:,j)=inv(A)*b(:,j);
end


n = zeros(N,11);
for j=1:11
for i =9:59
    n(i,j) = ni*exp(q*phi(i,j)/(k*T)); 
end
end


b2= zeros(N,11);
for j=1:11
b2(1,j)=0.33374-(j-1)/10; 
b2(N,j)=0.33374-(j-1)/10;
for i=interface1:interface2
    if (i==interface1) b2(i,j)= Delta_x*Delta_x*q*(Nacc+n(i,j))/eps0*0.5;
    elseif (i==interface2) b2(i,j)=  Delta_x*Delta_x*q*(Nacc+n(i,j))/eps0*0.5;
    else b2(i,j) = Delta_x*Delta_x*q*(Nacc+n(i,j))/eps0;
        end
end
phi2(:,j)=inv(A)*b2(:,j);
end
x=0.1*[0:(N-1)];
for j=1:11;
    for i=1:N
potential_diff(i,j)=phi(i,j)-phi2(i,j);
    end
end
plot(x,potential_diff);
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V');
