clear;

h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K

Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; % Lenghs, m
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0
nmax = 50;
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);

Ef = [-0.3:0.05:0.1];  % -0.3eV ~ 0.1eV

subbandNumber = 0;
totalNumber = zeros(length(Ef),1);


for i = 1:length(Ef)    
    for n=1:nmax
        Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;
        subbandNumber = coef*log(1+exp(-(Ez-Ef(i)*q)/(k_B*T)));
        totalNumber(i,1) = totalNumber(i,1) + subbandNumber;
    end    
end

semilogy(Ef,totalNumber/(Lx*Ly*1e4),'-o')  %/cm^2
xlabel('Fermi level (eV)')
ylabel('Integrated electron density (cm^{-2})')