clear;

q = 1.602192e-19; % Elementary charge, C
% eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
% eps_si = 11.7;
% eps_ox = 3.9;
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
ni = 1e16; % 1e10 /cm^3
Nd = zeros(9,1);
Na = zeros(9,1);

phi_Nd = zeros(9,1);
phi_analytic_Nd = zeros(9,1);
jaco1 = zeros(9,1);
res1  = zeros(9,1);
update1 = zeros(9,1);

for i = 1:9
    Nd(i,1) = 10^(15+i); %1e10 ~ 1e18 /cm^3
    if i>1
        phi_Nd(i,1)=phi_Nd(i-1,1);
    end
    for newton = 1:20;
        jaco1(i,1) = ni/thermal*exp(phi_Nd(i,1)/thermal)+ni/thermal*exp(-phi_Nd(i,1)/thermal);
        res1(i,1) = ni*exp(phi_Nd(i)/thermal)-ni*exp(-phi_Nd(i,1)/thermal)-Nd(i,1);
        update1(i,1) = -res1(i,1)/jaco1(i,1);
        phi_Nd(i,1) = phi_Nd(i,1) + update1(i,1);
        
        if abs(update1(i,1)) < 1e-10
            break;
        end
    end
    phi_analytic_Nd(i,1) = thermal*asinh(Nd(i,1)/2/ni);    
end

figure(1)
semilogx(Nd/1e6,phi_Nd,'o-'); hold on;
semilogx(Nd/1e6,phi_analytic_Nd,'rx');
xlabel('Nd (cm^{-3})')
ylabel('Electrostatic Potential (V)')
legend('Newton method', 'Analytic')

phi_Na = zeros(9,1);
phi_analytic_Na = zeros(9,1);
jaco2 = zeros(9,1);
res2  = zeros(9,1);
update2 = zeros(9,1);

for i = 1:9
    Na(i,1) = 10^(15+i); %1e10 ~ 1e18 /cm^3
    if i>1
        phi_Na(i,1)=phi_Na(i-1,1);
    end
    for newton = 1:20;
        jaco2(i,1) = ni/thermal*exp(phi_Na(i,1)/thermal)+ni/thermal*exp(-phi_Na(i,1)/thermal);
        res2(i,1) = ni*exp(phi_Na(i,1)/thermal)-ni*exp(-phi_Na(i,1)/thermal)+Na(i,1);
        update2(i,1) = -res2(i,1)/jaco2(i,1);
        phi_Na(i,1) = phi_Na(i,1) + update2(i,1);
        
        if abs(update2(i,1)) < 1e-10
            break;
        end
    end
    phi_analytic_Na(i,1) = thermal*asinh(-Na(i,1)/2/ni);
end

figure(2)
semilogx(-Na/1e6,phi_Na,'o-'); hold on;
semilogx(-Na/1e6,phi_analytic_Na,'rx');
xlabel('Na (cm^{-3})')
ylabel('Electrostatic Potential (V)')
legend('Newton method', 'Analytic')

figure(3)
error1 = zeros(9,1);
error1 = phi_Nd - phi_analytic_Nd;
semilogx(Nd/1e6,error1,'o-'); hold on;
xlabel('Nd (cm^{-3})')
ylabel('Difference (V)')

figure(4)
error2 = zeros(9,1);
error2 = phi_Na - phi_analytic_Na;
semilogx(-Na/1e6,error2,'o-'); hold on;
xlabel('Na (cm^{-3})')
ylabel('Difference (V)')
