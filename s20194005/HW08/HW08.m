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
coef = Delta_x*Delta_x*q/eps0;
thermal = k*T/q;
Jaco=sparse(N,N);
Jaco(1,1)=1;
Jaco(N,N)=1;
phi=zeros(N,11);
for i = 1:11
    phi(1,i) = 0.33374-(i-1)/10;
    phi(N,i) = 0.33374-(i-1)/10;
end
res= zeros(N,11);

for newton = 1:100
for j = 1:11
        for i = 2:66
            if  (i< interface1 || i> interface2)
                res(i,j) = eps_ox*phi(i+1,j)-2*eps_ox*phi(i,j)+eps_ox*phi(i-1,j);
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -2*eps_ox; 
                Jaco (i,i+1) = eps_ox;
            elseif  (i == interface1)
                res(i,j) = eps_si*phi(i+1,j)-(eps_si+eps_ox)*phi(i,j)+eps_ox*phi(i-1,j);
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -(eps_si+eps_ox); 
                Jaco(i,i+1) = eps_si;
            elseif  (i == interface2)
                res(i,j) = eps_ox*phi(i+1,j)-(eps_ox+eps_si)*phi(i,j)+eps_si*phi(i-1,j);
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -(eps_ox+eps_si); 
                Jaco(i,i+1) = eps_ox;
            else
                res(i,j) = eps_si*phi(i+1,j)-2*eps_si*phi(i,j)+eps_si*phi(i-1,j);
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -2*eps_si; 
                Jaco(i,i+1) = eps_si;
            end
        end
        
   for i = interface1:interface2
            if  (i == interface1)
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal))*0.5;
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal*0.5;
            elseif  (i==interface2)
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal))*0.5;
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal*0.5;
            else
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal));
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal;
            end
        end
        update(:,j) = Jaco\(-res(:,j));
        phi(:,j)=phi(:,j)+update(:,j);   
    end
end
density = zeros(N,11);
for j = 1:11
    for i = interface1:interface2
        density(i,j) = ni*exp(phi(i,j)/(thermal));
    end
end
density_sum=sum(density*0.1e-9,1);
x=0.1*[0:(N-1)];
plot(x,phi);
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1.0V');
semilogy(Applied_Voltage,density_sum*1e-4);
