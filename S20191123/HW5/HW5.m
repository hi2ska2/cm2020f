clear;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7;
eps_ox = 3.9;

Deltax = 0.1e-9; % 0.1 nm spacing
N = 67;
x = Deltax*transpose([0:N-1]); % real space, m
interface1 = 9; % At x=0.8 nm
interface2 = 59; % At x=5.8 nm
Nacc = [1e23 1e24 1e25]; % 1e17, 1e18, 1e19 /cm^3
coef = Deltax*Deltax*q/eps0;
A = zeros(N,N);
A(1,1) = 1;
A(N,N) = 1;

b = zeros(N,1);
X = zeros(N,3);
Analytic = zeros(N,3);

for k = 1:length(Nacc)
    for i = 2:N-1
        if(i < interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        elseif(i == interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*Nacc(k)*0.5;
        elseif(i == interface2)
            A(i,i-1) = eps_si;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_ox;
            b(i,1) = coef*Nacc(k)*0.5;
        elseif(i>interface2)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        else
            A(i,i-1) = eps_si;
            A(i,i) = -2*eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*Nacc(k);
        end
    end
    
    X(:,k) = A\b;  
    
    %Analytic
    for i = 1:interface1
        Analytic(i,k) = x(i)*(-3*0.8e-9*q*Nacc(k)*5e-9/2/(eps0*eps_si)/0.8e-9);
    end
    for i = interface2:N
        Analytic(i,k) = (x(i)-6.6e-9)*(3*0.8e-9*q*Nacc(k)*5e-9/2/(eps0*eps_si)/0.8e-9);
    end
    for i = interface1+1:interface2-1
        Analytic(i,k) = 1/2*q*Nacc(k)/(eps0*eps_si)*(x(i)-0.8e-9-2.5e-9)^2-3*0.8e-9*q*Nacc(k)*5e-9/2/(eps0*eps_si)-q*Nacc(k)/2/(eps0*eps_si)*(5e-9/2)^2;
    end
    
    
 
end

plot(x*1e9,X); hold on;
plot(x*1e9,Analytic,'o');
legend('1e17/cm^3', '1e18/cm^3', '1e19/cm^3','1e17/cm^3', '1e18/cm^3', '1e19/cm^3')
xlabel('Position (nm)')
ylabel('Potential (V)')
