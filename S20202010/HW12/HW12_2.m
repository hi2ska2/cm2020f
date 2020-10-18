%Assignment 12
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
%Nacc = 0; %undoped
Nacc = 1e24; %number density, m^-3, 10^18cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
Np1 = 1.0e16;
m0 = 9.1093837015e-31 ;
mxx = 0.19*m0 ; mzz = 0.19*m0; myy = 0.91*m0; %masses, m0
Lx = 100e-9; Ly = 100e-9; Lz=5e-9; %lengths, m
hbar = 1.054571800e-34; %Planck's constant/2pi, Js
coef2 = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Ltot = tox+tsi;



interface1 = 9; %at x =0.8nm
interface2 = 59; %at x =5.8nm
j1=interface1;
j2=interface2;
Deltax = (6.6e-9)/(N-1); %spacing
coef = Deltax*Deltax*q/eps0;


%%%%%%%%%%%%%%%%%%% Semiclassical electron density %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Numerical solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%% Boundary Condition %%%%%%%%%%%%%%%%%%%%%%%
    
    
phi = zeros(N,11);
phi(:,1) = 0.33374;
res = zeros (N,1);
Jaco = sparse (N,N);
res(1,1) = phi(1,1) - 0.33374;
res(N,1) = phi(N,1) - 0.33374;
GV = zeros(11,1);
elec = zeros(N,11);
int_elec = zeros(101,1); %integration of electron density
totalNumber = 0;
Nz = 51;
z = transpose([0:Nz-1])*Lz/(Nz-1);
elec_semi = zeros(Nz,1); % Electron density, /m^3
deltaz = Lz/(Nz-1);

for gate = 1:11
    if (gate ~=1)
        phi(:, gate)= phi(:, gate-1); %Using the previous solution for next gate Voltage
    end
    phi(1, gate) = 0.33374 + ((gate-1.0)/10.0);
    phi(N, gate) = phi(1, gate);
    Jaco(1,1) = 1.0;
    Jaco(N,N) = 1.0;
    res(1,1) = phi(1,gate) - 0.33374-((gate-1.0)/10.0);
    res(N,1) = res(1,1);
    GV(gate,1) = (gate-1.0)/10.0;
    
    

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

end
elec2 = zeros (51,11);
for gate = 1:11
    for ii = j1 : j2
        elec2(ii-j1+1,gate) =ni*exp(phi(ii,gate)/V_T);
    end
end
figure(1)
for gate = 1:11
    plot (z/1e-9, elec2(:,gate)/1e6, 'color', rand(1,3));
    hold on
     if gate == 11
        hold off
    end
    
end
xlabel('position x(nm)');
ylabel('Charge density (/cm^3)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');

%%%%%%%%%%%%%% quantum mechanical electron density %%%%%%%%%%%%%%
Vz = zeros(Nz,11);
Vz(:,:) = -q*phi(j1:j2,:)+q*0.56; %-q*phi +q*(Ec-Ei)
for gate = 1:11
    H = zeros(Nz-2,Nz-2);
    H(1,1) = -(-2)+Vz(2,gate)/(hbar^2/(2*mzz*deltaz^2));
    H(1,2) =  -1;
    for i = 2:Nz-3
        H(i,i-1) =  -1;
        H(i,i) = -(-2)+Vz(i+1,gate)/(hbar^2/(2*mzz*deltaz^2));
        H(i,i+1) =  -1;
    end
    H(Nz-2, Nz-2) = -(-2)+Vz(Nz-1,gate)/(hbar^2/(2*mzz*deltaz^2));
    H(Nz-2,Nz-3) =  -1;
    [V, D] = eig(H);
    D=diag(D)*(hbar^2/(2*mzz*deltaz^2)); %energy eigenvalue
    E=zeros(1,49);
    C=[E;V;E];
    C = C*sqrt(1/deltaz); %normalization


    n = zeros(Nz,11);
    for i = 1:Nz

        for j = 1:Nz-2
            n(i,gate) =n(i,gate)+ 2*C(i,j)*C(i,j)/(Lx*Ly)/(1+exp(D(j,1)/(k_B*T))); %n(z) 그래프 숫자 2 는 spin degeneracy 를 의미
        end
    end
    figure(2)
    plot (z/1e-9, n(:,gate)/1e6,'color', rand(1,3));
    hold on
    if gate == 11
        hold off
    end
end
xlabel('position x(nm)');
ylabel('Charge density (/cm^3)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');
    
    
    
    
    
    
    
 