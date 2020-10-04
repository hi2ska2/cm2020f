clear;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
Deltax = 0.1e-9; % 0.1 nm spacing
N = 67; % 6.6 nm thick
x = Deltax*transpose([0:N-1]); % real space, m
interface1 = 9; % At x=0.8 nm
interface2 = 59; % At x=5.8 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Nacc = 1e24; % 1e18 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltax*Deltax*q/eps0;




electron=zeros(N,1);
Vg=[0:0.02:1];
phi = zeros(N,length(Vg));
phi(:) = 0.33374;
int_electron = zeros(length(Vg),1);

for i = 1:length(Vg)   
    
    if i>1
        phi(:,i) = phi(:,i-1);
    end
        
    for newton = 1: 10
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
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,i)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,i)/thermal)/thermal*0.5;
            elseif (ii==interface2)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,i)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,i)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,i)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,i)/thermal)/thermal;
            end
        end
        update = Jaco \ (-res);                   
        phi(:,i) = phi(:,i) + update;
        
        error = max(abs(update));
        if error < 1e-10
            break;
        end
        
    end
    
    %%% integrate electron density %%%
    for ii=interface1:interface2
        if ii == interface1 || ii== interface2
            electron(ii,1) = ni*exp(phi(ii,i)/thermal)*0.5;
        else
            electron(ii,1) = ni*exp(phi(ii,i)/thermal);
        end
    end
    int_electron(i)=sum(electron)*Deltax;  %%/m^2
    
end
 
plot(Vg,int_electron/1e4,'o-'); hold on;  %%% linear scale
xlabel('Gate voltage (V)')
ylabel('Integrated electron density (cm^{-2})')

figure(2)
semilogy(Vg,int_electron/1e4,'o-'); hold on;   %%% log scale
xlabel('Gate voltage (V)')
ylabel('Integrated electron density (cm^{-2})')

figure(3)
plot(x*1e9,phi); hold on;
xlabel('Position (nm)')
ylabel('Electrostatic Potential (V)')