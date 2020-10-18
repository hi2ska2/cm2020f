clear;

h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
Lx = 100e-9; Ly = 100e-9; % Lenghs, m
Deltaz = 0.1e-9; % 0.1 nm spacing
N = 67; % 6.6 nm thick
z = Deltaz*transpose([0:N-1]); % real space, m
interface1 = 9; % At x=0.8 nm
interface2 = 59; % At x=5.8 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Nacc = 1e24; % 1e18 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
coef_Poi = Deltaz*Deltaz*q/eps0;
Ec_Ei = 0.56; % E_c-E_i, eV


Vg=[0:0.1:1];
electron_Poi=zeros(N,length(Vg));   % Electron density, /m^3
electron_Sch = zeros(N,length(Vg)); % Electron density, /m^3
V = zeros(N,length(Vg));            % Potential energy, J
phi = zeros(N,length(Vg));
phi(:) = 0.33374;

nSubband = 10;

for i = 1:length(Vg)
%     totalNumber = 0;
    
    %%% Non-linear Poisson equation %%%
    if i>1
        phi(:,i) = phi(:,i-1);
    end
    
    for newton = 1: 20
        res = zeros(N,1);
        Jaco = sparse(N,N);
        res(1,1) = phi(1,i) - 0.33374-Vg(i);
        Jaco(1,1) = 1.0;
        res(N,1) = phi(N,i) - 0.33374-Vg(i);
        Jaco(N,N) = 1.0;
        
        %%% Laplacian %%%
        for ii=2:N-1
            if (ii< interface1 || ii> interface2)
                res(ii,1) = eps_ox*phi(ii+1,i) - 2*eps_ox*phi(ii,i) + eps_ox*phi(ii-1,i);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco(ii,ii+1) = eps_ox;
            elseif (ii==interface1)
                res(ii,1) = eps_si*phi(ii+1,i) - (eps_si+eps_ox)*phi(ii,i) + eps_ox*phi(ii-1,i);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_si+eps_ox); Jaco(ii,ii+1) = eps_si;
            elseif (ii==interface2)
                res(ii,1) = eps_ox*phi(ii+1,i) - (eps_ox+eps_si)*phi(ii,i) + eps_si*phi(ii-1,i);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,i) - 2*eps_si*phi(ii,i) + eps_si*phi(ii-1,i);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
            end
        end
        
        %%% Charge %%%
        for ii=interface1:interface2
            if (ii==interface1)
                res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,i)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,i)/thermal)/thermal*0.5;
            elseif (ii==interface2)
                res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,i)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,i)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,i)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,i)/thermal)/thermal;
            end
        end
        update = Jaco \ (-res);
        phi(:,i) = phi(:,i) + update;
        
        error = max(abs(update));
        if error < 1e-10
            break;
        end
        
    end
    
    %%% electron density %%%
    for ii=interface1:interface2
        if ii == interface1 || ii== interface2
            electron_Poi(ii,i) = ni*exp(phi(ii,i)/thermal)/1e6; %/cm^3
        else
            electron_Poi(ii,i) = ni*exp(phi(ii,i)/thermal)/1e6; %/cm^3
        end
    end      
    
    %%% Schrhodinger eqaution %%%
    V(:,i) = q*Ec_Ei - q*phi(:,i); % Potential energy, J
    
    
    for iValley = 1:3
        mass = ones(3)*0.19;
        mass(iValley) = 0.91;
        coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mass(1)*mass(2))*m0/(hbar^2)*(k_B*T); %spin
        
        
        Nbulk = interface2-interface1-1; % Number of bulk silicon nodes
        Hamil = zeros(Nbulk,Nbulk);
        Hamil(1,1) = -2; Hamil(1,2) = 1;
        
        for ii=2:Nbulk-1
            Hamil(ii,ii+1) = 1;
            Hamil(ii,ii ) = -2;
            Hamil(ii,ii-1) = 1;
        end
        
        Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
        
        
        for ii=1:Nbulk
            Hamil(ii,ii) = Hamil(ii,ii) -2*mass(3)*m0*(Deltaz/hbar)^2*V(ii+interface1,i);
        end
        
        
        [eigenvectors,eigenvalues] = eig(Hamil);
        scaledEz = diag(eigenvalues)/(-2*mass(3)*m0*(Deltaz/hbar)^2); % Eigenenergy, J
        [sortedEz,sortedIndex] = sort(scaledEz);
        
        
        for n=1:nSubband
            Ez = sortedEz(n,1);
            wavefunction2 = eigenvectors(:,sortedIndex(n)).^2;
            normalization = sum(wavefunction2)*Deltaz;
            wavefunction2 = wavefunction2 / normalization;
            subbandNumber = 2*coef_Sch*log(1+exp(-Ez/(k_B*T))); % valley degeneracy
%             totalNumber = totalNumber + subbandNumber;
            electron_Sch(interface1+1:interface2-1,i) = electron_Sch(interface1+1:interface2-1,i) + 1/(Lx*Ly)*wavefunction2*subbandNumber;
        end
        
    end
end


electron_Sch = electron_Sch/1e6; %/cm^3

plot(z*1e9,electron_Poi,'r'); hold on;
plot(z*1e9,electron_Sch,'b'); hold on;
xlabel('z (nm)')
ylabel('electron density (cm^{-3})')