%Assignment 15_Long Channel(600nm)
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
un = 1500e-4; %electron mobility in silicon at 300K(1500cm^2/Vs)
m0 = 9.1093837015e-31 ;
mxx = 0.19*m0 ; myy = 0.19*m0; mzz = 0.91*m0; %masses, m0
Lx = 100e-9; Ly = 100e-9; Lz=6.6e-9; %lengths, m
hbar = 1.054571800e-34; %Planck's constant/2pi, Js
coef2 = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)/(hbar^2)*(k_B*T);




interface1 = (1/6)*(N-1); %at x =100 nm
interface2 = (5/6)*(N-1);
j1=interface1;
j2 = interface2;
Deltax = (600e-9)/(N-1); %spacing
coef = Deltax*Deltax*q/eps0;
x = Deltax * transpose([0:N-1]);
%%%%%%%%%%%%%%% Non_linear_Poisson's_equation %%%%%%%%%%%%%%%%%%%%%%
Ndop = 2e21*ones(N,1);
Ndop(1:j1,1) = 5e23;
Ndop(j2:N,1) = 5e23;
phi = zeros(N,1);
phi(:,1) = V_T*log(Ndop(:,1)/ni);
elec = zeros(N,1);
elec = ni*exp(phi/V_T);
res1=zeros(N,1);
Jaco1=sparse(N,N);


    
for nt = 1:10
    res1(1,1) = phi(1,1) - V_T*log(Ndop(1,1)/ni);
    Jaco1(1,1) = 1.0;


    for ii = 2:N-1
        res1(ii,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1));
        Jaco1(ii,ii+1) = eps_si;
        Jaco1(ii,ii) = -2*eps_si;
        Jaco1(ii,ii-1) = eps_si;
    end

    res1(N,1) = phi(N,1) - V_T*log(Ndop(1,1)/ni);
    Jaco1(N,N) = 1.0;

    for ii = 2:N-1
        
     
            res1(ii,1) = res1(ii,1) - coef*(ni*exp(phi(ii,1)/V_T) - Ndop(ii,1));
            Jaco1(ii,ii) = Jaco1(ii,ii) - coef*(ni*exp(phi(ii,1)/V_T)/V_T);
    end
        
    update1 = Jaco1\(-res1);
   phi = phi+update1;
end
figure(1)
elec1 = ni*exp(phi/V_T);
plot(x,elec1/1e6,'-');
hold on;




Ix = zeros(11,1);
Vx= zeros(11,1);

%%%% Scaling %%%%%
for bias = 1:11
    V_applied = 0.05*(bias-1);
    Vx(bias,1) = V_applied;

    for nt = 1:30
        res = zeros (2*N,1);
        Jaco = sparse(2*N,2*N);


        res(1,1) = phi(1,1) - V_T*log(Ndop(1,1)/ni);
        Jaco(1,1) = 1.0;


        for ii = 2:N-1
            res(2*ii-1,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1))+coef*(Ndop(ii,1)-elec(ii,1));
            Jaco(2*ii-1,2*ii+1) = eps_si;
            Jaco(2*ii-1,2*ii-1) = -2*eps_si;
            Jaco(2*ii-1,2*ii-3) = eps_si;
            Jaco(2*ii-1,2*ii) = -coef;
        end

        res(2*N-1,1) = phi(N,1) - V_T*log(Ndop(N,1)/ni)- V_applied;
        Jaco(2*N-1,2*N-1) = 1.0;


        for ii = 1:N-1 
            n_av = 0.5*(elec(ii+1,1)+elec(ii,1));


            dphidx = (phi(ii+1,1)-phi(ii,1))/Deltax;
            delecdx = (elec(ii+1,1) - elec(ii,1))/Deltax;



            Jn = n_av*dphidx - V_T*delecdx;


            res(2*ii,1) = res(2*ii,1) + Jn;
            Jaco(2*ii,2*ii+2) = Jaco(2*ii,2*ii+2)+0.5*dphidx-V_T/Deltax;
            Jaco(2*ii,2*ii)   = Jaco(2*ii,2*ii) +0.5*dphidx+V_T/Deltax;
            Jaco(2*ii, 2*ii+1) = Jaco(2*ii, 2*ii+1) + n_av/Deltax;
            Jaco(2*ii, 2*ii-1) = Jaco(2*ii, 2*ii-1) - n_av/Deltax;

            res(2*ii+2,1) = res(2*ii+2,1) - Jn;
            Jaco(2*ii+2,2*ii+2) = Jaco(2*ii+2,2*ii+2) - 0.5*dphidx + V_T/Deltax;
            Jaco(2*ii+2,2*ii)   = Jaco(2*ii+2,2*ii) - 0.5*dphidx - V_T/Deltax;
            Jaco(2*ii+2, 2*ii+1) = Jaco(2*ii+2, 2*ii+1) - n_av/Deltax;
            Jaco(2*ii+2, 2*ii-1) = Jaco(2*ii+2, 2*ii-1) + n_av/Deltax;


        end
        res(2,1) = elec(1,1) - Ndop(1,1);
        Jaco(2, :) = 0.0;
        Jaco(2,2) = 1.0;
        res(2*N,1) = elec(N,1) - Ndop(N,1);
        Jaco(2*N,:) = 0.0;
        Jaco(2*N,2*N) = 1.0;

        Cvector = zeros(2*N,1);
        Cvector(1:2:2*N-1,1) = V_T;
        Cvector(2:2:2*N,1) = max(abs(Ndop));
        Cmatrix = spdiags(Cvector, 0 , 2*N, 2*N);

        Jaco_scaled = Jaco*Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,2*N,2*N);
        Jaco_scaled = Rmatrix*Jaco_scaled;
        res_scaled = Rmatrix*res;
        update_scaled = Jaco_scaled\(-res_scaled);
        update = Cmatrix * update_scaled;

        phi = phi + update(1:2:2*N-1,1);
        elec = elec + update(2:2:2*N,1);
        norm(update(1:2:2*N-1),inf);


    end
    
    Ix(bias,1) = q*un*(elec(N,1)*(phi(N,1)-phi(N-1,1))/Deltax-V_T*(elec(N,1)-elec(N-1,1))/Deltax);


end
plot(x, elec/1e6, 'o');
xlabel('position(m)');
ylabel('carrier density(cm^-3)');

hold off;

figure(2)
plot(Vx,Ix/1e6);
xlabel('bias voltage(V)');
ylabel('Terminal current per unit area(\muA/\mum^2)');

