% Parameters
q=1.602192e-19;  % Elementary charge, C
h=6.626176e-34; hbar=h/(2*pi);   % Planck Constant, J s
m0=9.109543e-31; % Electron mass, kg
k_B=1.380662e-23; % Planck constant, J/K
T=300; % Temperature, K 

%Device informations
Lx=100e-9; Ly=100e-9; Lz=5e-9;
nmax=2;
mzz=0.91*m0; mxx=0.19*m0; myy=0.19*m0;  % Directional effective mass

elec=zeros(nmax,5);
dEF=0.05;
for iEF=1:9
    EF=(iEF-1)*dEF-0.3;
    for n=1:nmax
        Ezn=hbar^2*n^2*pi^2/2/mzz/Lz^2;
        elec(n,iEF)=Lx*Ly/(2*pi)*mxx/(hbar^2)*k_B*T*log(1+exp((-Ezn+EF*q)/k_B/T)); % electrons at each subbband
    end
end
eDensity=2*sum(elec)/(Lx*Ly*1e+4);  % total electron density (cm-2). Spin degeneracy is considered.

semilogy(-0.3:dEF:0.1,eDensity,'-o','linewidth',1.25);
ylabel('Electron density(/cm^2)');
xlabel('Fermi Level(eV)');
