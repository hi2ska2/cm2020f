clear;

eps1 = 3.9;
eps2 = 22;
eps0 = 8.854187817e-12; %F/m

N = 25; %2.4nm
interface = 6; %x=5nm
dx = 0.1e-9; %delta_x=0.1nm
l = [0:dx:2.4e-9];
A = zeros(N,N);
A(1,1) = 1;
A(N,N) = 1;

for i = 2:N-1
    if(i < interface)
        A(i,i-1) = eps1;
        A(i,i) = -2*eps1;
        A(i,i+1) = eps1;
    elseif(i == interface)
        A(i,i-1) = eps1;
        A(i,i) = -eps1 -eps2;
        A(i,i+1) = eps2;
    else
        A(i,i-1) = eps2;
        A(i,i) = -2*eps2;
        A(i,i+1) = eps2;
    end
end

b = zeros(N,1);
b(N,1) = 1; %phi(x=2.4nm)=1V

X = A\b;
plot(l*1e9,X,'DisplayName','numerical' );  hold on;


%analytic
X_analytical = zeros(length(X),1);

for i = 1:interface
    X_analytical(i) = 1/(2.4e-9*(5/24+19/24*3.9/22))*l(i);
    
end

for i = interface+1:length(X)
    X_analytical(i) = 3.9/22/(2.4e-9*(5/24+19/24*3.9/22))*l(i)+1-3.9/22/(2.4e-9*(5/24+19/24*3.9/22))*2.4e-9;
    
end
plot(l*1e9,X_analytical,'o','DisplayName','analytic');


xlabel('Position (nm)')
ylabel('Potential (V)')


%%% Capacitance %%%

D1 = eps1*eps0*(X(interface,1)-X(1,1))/0.5e-9; %eps1*delta_phi/0.5nm
D2 = eps2*eps0*(X(N,1)-X(interface,1))/1.9e-9; %eps2*delta_phi/1.9nm


C_numerical =D1/1*1e-4; %F/cm^2

C_1 = eps1/0.5e-9; C_2 = eps2/1.9e-9;
C_ananlytic = C_1*C_2/(C_1+C_2)*eps0*1e-4; %F/cm^2