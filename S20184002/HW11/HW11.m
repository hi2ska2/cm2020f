clear;
h = 6.626e-34;
hbar = h/(2*pi);
q = 1.6e-19;
m0 = 9.1e-31;
k_B = 1.38e-23;
T = 300;
Lx = 100e-9; Ly = 100e-9; Lz =5e-9;
mxx = 0.19; myy = 0.19; mzz = 0.91;
nmax = 10;
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Ef = -0.3*q:0.01*q:0.1*q;
TotalNumber = zeros(length(Ef),1);

for i =1:length(Ef)
    for n=1:nmax
    Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;
    SubbandNumber = coef*log(1+exp(-(Ez-Ef(1,i))/(k_B*T)));
    TotalNumber(i,1) = TotalNumber(i,1)+SubbandNumber;
    end
end

semilogy(Ef/q,TotalNumber/(Lx*Ly*1e4));
xlabel('Fermi Level (eV)') 
ylabel('Integrated Electron Density (cm^-2)')

