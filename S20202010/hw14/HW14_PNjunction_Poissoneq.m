%Assignment 14_fullarea
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide

N = 601; % 600 nm thick
k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.075e16; % 1.075e10/cm^3, intrinsic carrier density
T=300; %temp. 300K
%Nacc = 0; %undoped
Ndon = 1e23; %extrinsic doping density, m^-3, 10^17cm^-3
Nacc = 1e23; %number density, m^-3, 10^17cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)

m0 = 9.1093837015e-31 ;
mxx = 0.19*m0 ; myy = 0.19*m0; mzz = 0.91*m0; %masses, m0
Lx = 100e-9; Ly = 100e-9; Lz=6.6e-9; %lengths, m
hbar = 1.054571800e-34; %Planck's constant/2pi, Js
coef2 = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)/(hbar^2)*(k_B*T);




interface1 = 301; %at x =300 nm

j1=interface1;

Deltax = (600e-9)/(N-1); %spacing
coef = Deltax*Deltax*q/eps0;
x = Deltax * transpose([0:N-1]);
%%%%%%%%%%%%%%% Non_linear_Poisson's_equation %%%%%%%%%%%%%%%%%%%%%%
Ndop = zeros(N,1);
Ndop(:,1) = 1e23;
res = zeros (N,1);
Jaco = sparse(N,N);
phi = zeros(N,1);
phi(1:j1,1) = -V_T*log(Ndop(1:j1,1)/ni); %p-type
phi(j1+1:N,1) = V_T*log(Ndop(j1+1:N,1)/ni); %n-type

for nt = 1:10
    res(1,1) = phi(1,1) + V_T*log(Ndop(1,1)/ni);
    Jaco(1,1) = 1.0;


    for ii = 2:N-1
        res(ii,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1));
        Jaco(ii,ii+1) = eps_si;
        Jaco(ii,ii) = -2*eps_si;
        Jaco(ii,ii-1) = eps_si;
    end

    res(N,1) = phi(N,1) - V_T*log(Ndop(1,1)/ni);
    Jaco(N,N) = 1.0;

    for ii = 2:N-1
        if (ii < j1+1)
            res(ii,1) = res(ii,1) + coef*(ni*exp(-phi(ii,1)/V_T) - ni*exp(phi(ii,1)/V_T) - Ndop(ii,1));
            Jaco(ii,ii) = Jaco(ii,ii) - coef*(ni*exp(-phi(ii,1)/V_T)/V_T + ni*exp(phi(ii,1)/V_T)/V_T);
  
        else
            res(ii,1) = res(ii,1) + coef*(ni*exp(-phi(ii,1)/V_T) - ni*exp(phi(ii,1)/V_T) + Ndop(ii,1));
            Jaco(ii,ii) = Jaco(ii,ii) - coef*(ni*exp(-phi(ii,1)/V_T)/V_T + ni*exp(phi(ii,1)/V_T)/V_T);
        end
    end

   update = Jaco\(-res);
   phi = phi+update;
   
    
end
figure(1);
elec = zeros(N,1);
elec = ni*exp(phi/V_T);
plot(x,elec/1e6,'b');
hold on;
hole = zeros(N,1);
hole = ni*exp(-phi/V_T);
plot(x, hole/1e6, 'r');
hold on;



%%%%%%%%%%%%%%%%%%% continuity equation %%%%%%%%%%%%%%%%%%%%%%

res_elec = zeros(N,1);
Jaco_elec = sparse(N,N);
res_hole = zeros(N,1);
Jaco_hole = sparse(N,N);

%boundary condition


%construct Jaco and res

for ii = 1:N-1 
    n_av = 0.5*(elec(ii+1,1)+elec(ii,1));
    h_av = 0.5*(hole(ii+1,1)+hole(ii,1));
    
    dphidx = (phi(ii+1,1)-phi(ii,1))/Deltax;
    delecdx = (elec(ii+1,1) - elec(ii,1))/Deltax;
    dholedx = (hole(ii+1,1)- hole(ii,1))/Deltax;
    
    
    Jn = n_av*dphidx - V_T*delecdx;
    Jp = h_av*dphidx +V_T*dholedx;
    
    res_elec(ii,1) = res_elec(ii,1) + Jn;
    Jaco_elec(ii,ii+1) = Jaco_elec(ii,ii+1)+0.5*dphidx-V_T/Deltax;
    Jaco_elec(ii,ii)   = Jaco_elec(ii,ii) +0.5*dphidx+V_T/Deltax;
    
    res_hole(ii,1) = res_hole(ii,1) + Jp;
    Jaco_hole(ii,ii+1) = Jaco_hole(ii,ii+1)+0.5*dphidx+V_T/Deltax;
    Jaco_hole(ii,ii)   = Jaco_hole(ii,ii) +0.5*dphidx-V_T/Deltax;
    
    res_elec(ii+1,1) = res_elec(ii+1,1) - Jn;
    Jaco_elec(ii+1,ii+1) = Jaco_elec(ii+1,ii+1)-0.5*dphidx+V_T/Deltax;
    Jaco_elec(ii+1,ii)   = Jaco_elec(ii+1,ii) -0.5*dphidx-V_T/Deltax;
    
    res_hole(ii+1,1) = res_hole(ii+1,1) - Jp;
    Jaco_hole(ii+1,ii+1) = Jaco_hole(ii+1,ii+1)-0.5*dphidx-V_T/Deltax;
    Jaco_hole(ii+1,ii)   = Jaco_hole(ii+1,ii) -0.5*dphidx+V_T/Deltax;
    

end
res_elec(1,1) = elec(1,1);
res_hole(1,1) = hole(1,1) - Ndop(1,1); %아리까리
Jaco_elec(1,:) = 0.0;
Jaco_elec(1,1) = 1.0;
Jaco_hole(1,:) = 0.0;
Jaco_hole(1,1) = 1.0;
res_elec(N,1) = elec(N,1) - Ndop(N,1);
res_hole(N,1) = hole(N,1);
Jaco_elec(N,:) = 0.0;
Jaco_elec(N,N) = 1.0;
Jaco_hole(N,:) = 0.0;
Jaco_hole(N,N) = 1.0;

update_elec = Jaco_elec\(-res_elec);
elec = elec + update_elec;
update_hole = Jaco_hole\(-res_hole);
hole = hole + update_hole;

plot (x, elec/1e6,'o','color', 'b');
hold on;
plot (x, hole/1e6,'o','color', 'r');
hold off;

xlabel('position(m)');
ylabel('carrier density(cm^-3)');
legend('Psolve electron','Psolve hole','Csolve electron','Csolve hole', 'Location','Bestoutside');


