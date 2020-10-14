%Assignment 11
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
m0 = 9.1093837015e-31 ; %electron rest mass, kg
T=300; %temp. 300K
Nacc = 1e26; %number density, m^-3, 10^20cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
mxx = 0.19 ; myy = 0.19; mzz = 0.91; %masses, m0
Lx = 100e-9; Ly = 100e-9; Lz=5e-9; %lengths, m
hbar = 1.054571800e-34; %Planck's constant/2pi, Js
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);

V_F = -0.3; %Fermi voltage, V
V=zeros(101,1);
nmax = 50; % 충분히 큰 수
totalNumber = zeros(101,1);
for z = 1:101;
    V(z,1)=V_F;
    
    for n = 1:nmax
        Ez = hbar^2/(2*mzz*m0)*(pi*n/Lz)^2;
        subbandNumber = coef*log(1+exp((-Ez+q*V_F)/(k_B*T)));
        totalNumber(z,1) = totalNumber(z,1) + subbandNumber;
    end
    V_F = V_F+0.4/100;
end

int_total = totalNumber/(Lx*Ly);
semilogy(V(:,1),int_total(:,1)*1e-4, 'o-');
xlabel('Fermi Energy (eV)');
ylabel('Integrated electron density (cm^-2)');