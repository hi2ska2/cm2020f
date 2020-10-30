q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
Deltax = 1e-9; % 1 nm spacing
N = 501; % 500-nm-long structure
x = Deltax*transpose([0:N-1]); % real space, m
interface = 251; % At x=250 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Ndon = 1e23*ones(N,1); % 1e17 /cm^3
% Ndon(1:interface,1) = -1e23; % P
ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltax*Deltax*q/eps0;


phi = zeros(N,1);
phi(1:interface,1) = -thermal*log(Ndon(1:interface,1)/ni);
phi(interface+1:N,1) = thermal*log(Ndon(interface+1:N,1)/ni);

for newton=1:10
    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,1) + thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1.0;
    for ii=2:N-1
        res(ii,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1));
        Jaco(ii,ii-1) = eps_si;
        Jaco(ii,ii ) = -2*eps_si;
        Jaco(ii,ii+1) = eps_si;
    end
    res(N,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(N,N) = 1.0;
    
    for ii=2:interface
        res(ii,1) = res(ii,1) - coef*(Ndon(ii,1)+ni*exp(phi(ii,1)/thermal)-ni*exp(-phi(ii,1)/thermal));
        Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,1)/thermal)/thermal - coef*ni*exp(-phi(ii,1)/thermal)/thermal;
    end
    
    for ii=interface+1:N-1
        res(ii,1) = res(ii,1) - coef*(-Ndon(ii,1)+ni*exp(phi(ii,1)/thermal)-ni*exp(-phi(ii,1)/thermal));
        Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,1)/thermal)/thermal - coef*ni*exp(-phi(ii,1)/thermal)/thermal;
    end
    update = Jaco \ (-res);
    phi = phi + update;
    norm(update,inf)
end
% plot(x,phi);

elec = zeros(N,1);
hole = zeros(N,1);
elec = ni*exp(phi/thermal);
hole = ni*exp(-phi/thermal);
figure(1)
semilogy(x*1e9,elec/1e6,'r'); hold on;
figure(2)
semilogy(x*1e9,hole/1e6,'r');hold on;



%%% continuity equation, electron %%%
res_elec = zeros(N,1);
Jaco_elec = sparse(N,N);

res_elec(1,1) = elec(1,1) - ni^2/Ndon(1,1);
Jaco_elec(1,:) = 0.0;
Jaco_elec(1,1) = 1.0;

for ii=2:N-1 
    res_elec(ii,1) = [0.5*(elec(ii+1,1)+elec(ii,1)) * (phi(ii+1,1)-phi(ii,1))/Deltax - thermal * (elec(ii+1,1)-elec(ii,1))/Deltax] - [0.5*(elec(ii,1)+elec(ii-1,1)) * (phi(ii,1)-phi(ii-1,1))/Deltax - thermal * (elec(ii,1)-elec(ii-1,1))/Deltax];
    Jaco_elec(ii,ii+1) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax - thermal / Deltax;
    Jaco_elec(ii,ii ) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax + thermal / Deltax - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax + thermal / Deltax;
    Jaco_elec(ii,ii-1) = - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax - thermal / Deltax;    
end


res_elec(N,1) = elec(N,1) - Ndon(N,1);
Jaco_elec(N,:) = 0.0;
Jaco_elec(N,N) = 1.0;

update_elec = Jaco_elec \ (-res_elec);
elec = elec + update_elec;




%%% continuity equation, hole %%%
res_hole = zeros(N,1);
Jaco_hole = sparse(N,N);

res_hole(1,1) = hole(1,1) - Ndon(1,1);
Jaco_hole(1,:) = 0.0;
Jaco_hole(1,1) = 1.0;

for ii=2:N-1 
    res_hole(ii,1) = [0.5*(hole(ii+1,1)+hole(ii,1)) * (phi(ii+1,1)-phi(ii,1))/Deltax + thermal * (hole(ii+1,1)-hole(ii,1))/Deltax] - [0.5*(hole(ii,1)+hole(ii-1,1)) * (phi(ii,1)-phi(ii-1,1))/Deltax + thermal * (hole(ii,1)-hole(ii-1,1))/Deltax];
    Jaco_hole(ii,ii+1) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax + thermal / Deltax;
    Jaco_hole(ii,ii ) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax - thermal / Deltax - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax - thermal / Deltax;
    Jaco_hole(ii,ii-1) = - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax + thermal / Deltax;    
end


res_hole(N,1) = hole(N,1) - ni^2/Ndon(N,1);
Jaco_hole(N,:) = 0.0;
Jaco_hole(N,N) = 1.0;

update_hole = Jaco_hole \ (-res_hole);
hole = hole + update_hole;

figure(1)
semilogy(x*1e9,elec/1e6,'o'); hold on;
xlabel('Position (nm)')
ylabel('electron density (cm^{-3})')

figure(2)
semilogy(x*1e9,hole/1e6,'o');hold on;
xlabel('Position (nm)')
ylabel('hole density (cm^{-3})')