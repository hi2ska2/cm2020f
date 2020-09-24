%Assignment 7, Nonlinear
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, m
tsi = 5e-9; %silicon layer thickness, m

k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.0e16; % 1.0e10/cm^3, intrinsic carrier density
T=300; %temp. 300K
Nacc = 1e24; %number density, m^-3, 10^18cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
Np1 = 1.0e16;
%%%%%%%%%%%%%%%%%%%%% Numerical solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Ver. 1 Newton method %%%%%%%%%

soln1 = zeros(100,1);
soln1(1,1) = 1;
x = zeros(100,1);
x(100,1) = 100;
solnA2 = zeros(100,1);
solnA2(100,1) = V_T*asinh(Np1/(2*ni));
for nt = 1:(100-1)
    
    res1 = Np1 + ni*exp(-(soln1(nt,1)/V_T)) - ni*exp(soln1(nt,1)/V_T);
    Jaco1 = -(ni/V_T)*exp(-soln1(nt,1)/V_T) - (ni/V_T)*exp(soln1(nt,1)/V_T);
    update = Jaco1\(-res1);
    soln1(nt+1,1) = soln1(nt,1) + update;
    x(nt,1) = nt;
    solnA2(nt,1) = V_T*asinh(Np1/(2*ni));
end
figure(1)
plot(x, soln1,'o', x, solnA2, '-');
xlabel('iteration (#)'); 
ylabel('Solution');
legend('Numerical solution','Analytical solution','Location','bestoutside');
%%%%%%%%%%% Ver. 2 Newton method with various Np %%%%%%%%


soln2 = zeros(9,1);
solnA = zeros(9,1);
Np2 = 1.0e16;
x2 = zeros(9,1);
for ntt= 1:9
    soln2(ntt,1)= 1;
    solnA(ntt,1) = V_T*asinh(Np2/(2*ni));
    for nt = 1:50
        
    
    res2 = Np2 + ni*exp(-(soln2(ntt,1))/V_T) - ni*exp(soln2(ntt,1)/V_T);
    Jaco2 = -(ni/V_T)*exp(-soln2(ntt,1)/V_T) - (ni/V_T)*exp(soln2(ntt,1)/V_T);
    update2 = Jaco2\(-res2);
    soln2(ntt,1) = soln2(ntt,1) + update2;
    
    
    end
    Np2 = Np2*10;
    err= (solnA(ntt,1)-soln2(ntt,1))*100/solnA(ntt,1);
    x2(ntt,1) = 9 + ntt;
end

figure(2)
plot(x2, soln2, 'o', x2, solnA, ':');
xlabel('N^+ (10^x /cm^3)'); 
ylabel('Solution');
legend('Numerical solution','Analytical solution','Location','bestoutside');

figure(3)
plot(x2, abs(err));
xlabel('N'); 
ylabel('error(%)');

