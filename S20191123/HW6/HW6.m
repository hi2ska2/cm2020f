clear;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7;
eps_ox = 3.9;

k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
ni = 1.075e16; % 1.075e10 /cm^3

Deltax = 0.1e-9; % 0.1 nm spacing
N = 67;
x = Deltax*transpose([0:N-1]); % real space, m
interface1 = 9; % At x=0.8 nm
interface2 = 59; % At x=5.8 nm
Nacc = 1e24; % 1e18cm^3
coef = Deltax*Deltax*q/eps0;
A = zeros(N,N);
A(1,1) = 1;
A(N,N) = 1;

b = zeros(N,1);

V = [0:0.1:1];
X_initial = zeros(N,length(V));
X_update = zeros(N,length(V));
elec_initial = zeros(N,length(V));
elec_update = zeros(N,length(V));

Difference = zeros(N,length(V));
Difference_max = zeros(length(V),1);

for v = 1:length(V);
    
    b(1,1) = 0.33374+V(v);  %workfunction = 4.3eV.
    b(N,1) = 0.33374+V(v);
    
    for i = 2:N-1
        if(i < interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        elseif(i == interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*Nacc*0.5;
        elseif(i == interface2)
            A(i,i-1) = eps_si;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_ox;
            b(i,1) = coef*Nacc*0.5;
        elseif(i>interface2)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        else
            A(i,i-1) = eps_si;
            A(i,i) = -2*eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*Nacc;
        end
    end
    
    X_initial(:,v) = A\b;
    
    
    for i = interface1 : interface2
        elec_initial(i,v) = ni*exp(X_initial(i,v)/thermal);        
    end
    
    
    %%% update %%%
    b(1,1) = 0.33374+V(v);
    b(N,1) = 0.33374+V(v);
    
    for i = 2:N-1
        if(i < interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        elseif(i == interface1)
            A(i,i-1) = eps_ox;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*(Nacc+elec_initial(i,v))*0.5;
        elseif(i == interface2)
            A(i,i-1) = eps_si;
            A(i,i) = -eps_ox -eps_si;
            A(i,i+1) = eps_ox;
            b(i,1) = coef*(Nacc+elec_initial(i,v))*0.5;
        elseif(i>interface2)
            A(i,i-1) = eps_ox;
            A(i,i) = -2*eps_ox;
            A(i,i+1) = eps_ox;
        else
            A(i,i-1) = eps_si;
            A(i,i) = -2*eps_si;
            A(i,i+1) = eps_si;
            b(i,1) = coef*(Nacc+elec_initial(i,v));
        end
    end
    
    X_update(:,v) = A\b;
    
    for i = interface1 : interface2
        elec_update(i,v) = ni*exp(X_update(i,v)/thermal);        
    end
    
    Difference(:,v) =  X_update(:,v) - X_initial(:,v);
    Difference_max(v,1) = max(abs(Difference(:,v)));
end
figure(1)
plot(x*1e9,elec_initial/1e6,'o'); hold on;
plot(x*1e9,elec_update/1e6); 
xlabel('Position (nm)')
ylabel('electron density (cm^{-3})')
legend('Vg=0V', 'Vg=0.1V', 'Vg=0.2V', 'Vg=0.3V', 'Vg=0.4V', 'Vg=0.5V', 'Vg=0.6V', 'Vg=0.7V', 'Vg=0.8V', 'Vg=0.9V', 'Vg=1V')

figure(2)
plot(x*1e9,X_initial,'o'); hold on;
plot(x*1e9,X_update); 
xlabel('Position (nm)')
ylabel('Electrostatic Potential (V)')
legend('Vg=0V', 'Vg=0.1V', 'Vg=0.2V', 'Vg=0.3V', 'Vg=0.4V', 'Vg=0.5V', 'Vg=0.6V', 'Vg=0.7V', 'Vg=0.8V', 'Vg=0.9V', 'Vg=1V')


figure(3)
plot(V,Difference_max)
xlabel('Gate voltage (V)')
ylabel('Difference (V)')

