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
% electron_Poi=zeros(N,length(Vg));   % Electron density, /m^3
electron_Sch = zeros(N,length(Vg)); % Electron density, /m^3
V = zeros(N,length(Vg));            % Potential energy, J
phi = zeros(N,length(Vg));
phi(:,1) = 0.33374;

nSubband = 10;

%%% Non-linear Poisson equation %%%
for newton = 1: 20
    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,1) - 0.33374;
    Jaco(1,1) = 1.0;
    res(N,1) = phi(N,1) - 0.33374;
    Jaco(N,N) = 1.0;
    
    %%% Laplacian %%%
    for ii=2:N-1
        if (ii< interface1 || ii> interface2)
            res(ii,1) = eps_ox*phi(ii+1,1) - 2*eps_ox*phi(ii,1) + eps_ox*phi(ii-1,1);
            Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco(ii,ii+1) = eps_ox;
        elseif (ii==interface1)
            res(ii,1) = eps_si*phi(ii+1,1) - (eps_si+eps_ox)*phi(ii,1) + eps_ox*phi(ii-1,1);
            Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_si+eps_ox); Jaco(ii,ii+1) = eps_si;
        elseif (ii==interface2)
            res(ii,1) = eps_ox*phi(ii+1,1) - (eps_ox+eps_si)*phi(ii,1) + eps_si*phi(ii-1,1);
            Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
        else
            res(ii,1) = eps_si*phi(ii+1,1) - 2*eps_si*phi(ii,1) + eps_si*phi(ii-1,1);
            Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
        end
    end
    
    %%% Charge %%%
    for ii=interface1:interface2
        if (ii==interface1)
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
        elseif (ii==interface2)
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
        else
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal));
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal;
        end
    end
    update = Jaco \ (-res);
    phi(:,1) = phi(:,1) + update;
    
    error = max(abs(update));
    if error < 1e-10
        break;
    end
    
%         %%% Poisson electron density %%%
%     for ii=interface1:interface2
%         if ii == interface1 || ii== interface2
%             electron_Poi(ii,1) = ni*exp(phi(ii,1)/thermal)/1e6; %/cm^3
%         else
%             electron_Poi(ii,1) = ni*exp(phi(ii,1)/thermal)/1e6; %/cm^3
%         end
%     end  
    
end


Schrhodinger_loop = 40;
Update_MAX = zeros(Schrhodinger_loop,length(Vg));
for i = 1:length(Vg)
    
    if i>1
        phi(:,i) = phi(:,i-1);
    end
    
    %%% Schrhodinger eqaution %%%
    for iNewton = 1:Schrhodinger_loop
        
        V(:,i) = q*Ec_Ei - q*phi(:,i); % Potential energy, J
        electron_Sch(:,i) = 0; % Electron density, /m^3
        %     totalNumber = 0;
        
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
            
            %%% Calculation of the electron density %%%
            
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
        
        %%% Poisson equation %%%
        
        for poisson_newton = 1:20
            
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
                    res(ii,1) = res(ii,1) - coef_Poi*(Nacc)*0.5;
                    Jaco(ii,ii) = Jaco(ii,ii);
                elseif (ii==interface2)
                    res(ii,1) = res(ii,1)- coef_Poi*(Nacc)*0.5;
                    Jaco(ii,ii) = Jaco(ii,ii);
                else
                    res(ii,1) = res(ii,1) - coef_Poi*(Nacc+electron_Sch(ii,i));
                    Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*electron_Sch(ii,i)/thermal;
                end
            end
            update = Jaco \ (-res);
            phi(:,i) = phi(:,i) + update;
            
            for ii=interface1+1:interface2-1
                electron_Sch(ii,i) = electron_Sch(ii,i).*exp(update(ii,1)/thermal);
            end
            
            
            if poisson_newton == 1
                Update_MAX(iNewton,i) = max(abs(update));
            end
            
            error = max(abs(update));
            if error < 1e-10
                break;
            end
        end
    end   
end


electron_Sch = electron_Sch/1e6; %/cm^3

plot(z*1e9,electron_Sch,'r'); hold on;
xlabel('z (nm)')
ylabel('electron density (cm^{-3})')