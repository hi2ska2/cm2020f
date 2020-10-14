hbar = 1.054571e-34;
q = 1.602192e-19;
m0 = 9.109534e-31;
T = 300;
k_B = 1.380662e-23;
Lx = 100e-9; 
Ly = 100e-9; 
Lz = 5e-9;
nmax = 50;
mxx=0.19;
myy=0.19;
mzz=0.91;
coef=2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
totalNumber=zeros(41,1);
for i=1:41
    Ef=0.01*i-0.31;
    for n=1:nmax
       Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2; 
       subbandNumber = coef*log(1+exp((-Ez+q*Ef)/(k_B*T)));
       totalNumber(i,1)= totalNumber(i,1)+ subbandNumber;
    end
end

semilogy(-0.3:0.01:0.1,totalNumber/(Lx*Ly)*1e-4, '-o');
xlabel('Fermi energy (eV)')
ylabel('Electron density (cm^-2)')
