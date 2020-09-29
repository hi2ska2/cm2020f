%Assignment 8
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, m
tsi = 5e-9; %silicon layer thickness, m
N = 67; % 6.6 nm thick
k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.0e16; % 1.0e10/cm^3, intrinsic carrier density
T=300; %temp. 300K
Nacc = 1e24; %number density, m^-3, 10^18cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
Np1 = 1.0e16;

interface1 = 9; %at x =0.8nm
interface2 = 59; %at x =5.8nm
j1=interface1;
j2=interface2;
Deltax = (6.6e-9)/(N-1); %spacing
coef = Deltax*Deltax*q/eps0;
%%%%%%%%%%%%%%%%%%%%% Numerical solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%% Boundary Condition %%%%%%%%%%%%%%%%%%%%%%%
    
    
phi = zeros(N,101);
phi(:,1) = 0.33374;
res = zeros (N,1);
Jaco = sparse (N,N);
res(1,1) = phi(1,1) - 0.33374;
res(N,1) = phi(N,1) - 0.33374;
GV = zeros(101,1);
elec = zeros(N,101);
int_elec = zeros(101,1); %integration of electron density

for gate = 1:101
    if (gate ~=1)
        phi(:, gate)= phi(:, gate-1); %Using the previous solution for next gate Voltage
    end
    phi(1, gate) = 0.33374 - ((gate-1.0)/100.0);
    phi(N, gate) = phi(N, gate);
    Jaco(1,1) = 1.0;
    Jaco(N,N) = 1.0;
    res(1,1) = phi(1,gate) - 0.33374+((gate-1.0)/100.0);
    res(N,1) = res(1,1);
    GV(gate,1) = (gate-1.0)/100.0;
    
    

    for nt = 1:10


    for ii = 2 : N-1
%%%%%%%%%%%%%%%%Laplacian part%%%%%%%%%%%%%%%%%%%%%
            if (ii < j1 || ii>j2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - 2*eps_ox*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco(ii,ii+1) = eps_ox;
            elseif (ii==j1)
                res(ii,1) = eps_si*phi(ii+1,gate) - (eps_si+eps_ox)*phi(ii,gate) + eps_ox*phi(ii-1,gate) ;
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_si;
            elseif (ii==j2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - (eps_si+eps_ox)*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,gate) - 2*eps_si*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
            end
        end
%%%%%%%%%%%%%%%%%charge part%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii = j1:j2
            if(ii==j1)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/V_T))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/V_T)/V_T*0.5;
            elseif (ii==j2)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/V_T))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/V_T)/V_T*0.5;
            else
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/V_T));
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/V_T)/V_T;
            end
        end

        update = Jaco\(-res);
        phi(:,gate) = phi(:,gate) + update; 
    end 
 %%%%%%%%%%%%%%%%integration part%%%%%%%%%%%%%%%%%%%
   for ii = j1:j2
       elec(ii,gate) = ni*exp(phi(ii,gate)/V_T);
       if (ii==j1 || ii==j2)
           int_elec(gate,1) = int_elec(gate,1) + Deltax/2*elec(ii,gate);
       else
           int_elec(gate,1) = int_elec(gate,1) + Deltax/2*elec(ii,gate);
       end
   end
end

figure(1)
semilogy(GV(:,1),int_elec(:,1)*1e-4, 'o-');
xlabel ('Gate voltage (V)');
ylabel ('Integrated electron density (cm^-^2)');


