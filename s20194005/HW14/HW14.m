clear;
q = 1.602192e-19; 
eps0 = 8.854187e-12; 
Deltax = 1e-9;
m0 = 9.109383e-31 ;
N = 501; 
T=300.0;
hbar = 1.054571e-34;
k = 1.380662e-23;
ni = 1.075e16;
eps_si = 11.7; 
eps_ox = 3.9; 
Lx=100e-9;
Ly=100e-9;
x = Deltax * transpose([0:N-1]);
thermal = k*T/q;
coef = Deltax*Deltax*q/eps0;
Nd= zeros(N,1);
Nd(:,1) = 1e23;
interface = 251;
res = zeros (N,1);
Jaco = sparse (N,N);
phi = zeros (N,1);
phi(1:interface,1)=-thermal*log(Nd(1:interface,1)/ni);
phi(interface+1:N,1)=thermal*log(Nd(interface+1:N,1)/ni);
for newton=1:10

    res(1,1)= phi(1,1) + thermal*log(Nd(1,1)/ni);
    Jaco(1,1) =1;
    for j=2:N-1
        res(j,1)= eps_si*(phi(j+1,1)-2*phi(j,1)+phi(j-1,1));
        Jaco(j,j+1) = eps_si;
        Jaco(j,j) = -2*eps_si;
        Jaco(j,j-1) = eps_si;
    end
    res(N,1) = phi(N,1) - thermal*log(Nd(1,1)/ni);
    Jaco(N,N) = 1;
    
    for j=2:interface
    res(N,1) = res(j,1)-coef*(Nd(j,1)+ni*exp(phi(j,1)/thermal)-ni*exp(-phi(j,1)/thermal));
    Jaco(j,j)= Jaco(j,j)-coef*ni*exp(phi(j,1)/thermal)/thermal- coef*ni*exp(-phi(j,1)/thermal)/thermal;
    end
    for j=interface+1:N-1
    res(N,1) = res(j,1)-coef*(-Nd(j,1)+ni*exp(phi(j,1)/thermal)-ni*exp(-phi(j,1)/thermal));
    Jaco(j,j)= Jaco(j,j)-coef*ni*exp(phi(j,1)/thermal)/thermal- coef*ni*exp(-phi(j,1)/thermal)/thermal;
    end
    update = Jaco\(-res);
    phi = phi + update;
end 
figure(1);
elec = zeros(N,1);
elec = ni*exp(phi/thermal);
hole = zeros(N,1);
hole = ni*exp(-phi/thermal);
plot(x,elec/1e6,'color', 'r');
hold on
plot(x, hole/1e6,'color', 'b');
hold on;

%continuity equation

res_elec = zeros(N,1);
res_hole = zeros(N,1);
Jaco_elec = sparse(N,N);
Jaco_hole = sparse(N,N);
%for electron
for j=1:N-1
    n_av=0.5*(elec(j+1,1)+elec(j,1));
    dphidx= (phi(j+1,1)-phi(j,1))/Deltax;
    delecdx=(elec(j+1,1)-elec(j,1))/Deltax;
    Jn= n_av * dphidx- thermal *delecdx;
    res_elec(j,1)=res_elec(j,1)+ Jn;
    Jaco_elec(j,j+1)= Jaco_elec(j,j+1)+0.5*dphidx+thermal/Deltax;
    Jaco_elec(j,j)= Jaco_elec(j,j)+0.5*dphidx+thermal/Deltax;
    res_elec(j+1,1)= res_elec(j+1,1)-Jn;
    Jaco_elec(j+1,j+1)=Jaco_elec(j+1,j+1)-0.5*dphidx+thermal/Deltax;
    Jaco_elec(j+1,j)= Jaco_elec(j+1,j)-0.5*dphidx-thermal/Deltax;
end
 
res_elec(1,1) = elec(1,1);
Jaco_elec(1,:) = 0.0;
Jaco_elec(1,1) = 1.0;    
res_elec(N,1) = elec(N,1) - Nd(N,1);
Jaco_elec(N,:) = 0.0;
Jaco_elec(N,N) = 1.0;
update_elec = Jaco_elec\(-res_elec);
elec= elec + update_elec;

%for hole
for j=1:N-1
    p_av=0.5*(hole(j+1,1)+hole(j,1));
    dphidx= (phi(j+1,1)-phi(j,1))/Deltax;
    dholedx=(hole(j+1,1)-hole(j,1))/Deltax;
    Jp= p_av * dphidx- thermal *dholedx;
    res_hole(j,1)=res_hole(j,1)+ Jp;
    Jaco_hole(j,j+1)= Jaco_hole(j,j+1)+0.5*dphidx+thermal/Deltax;
    Jaco_hole(j,j)= Jaco_hole(j,j)+0.5*dphidx+thermal/Deltax;
    res_hole(j+1,1)= res_hole(j+1,1)-Jn;
    Jaco_hole(j+1,j+1)=Jaco_hole(j+1,j+1)-0.5*dphidx+thermal/Deltax;
    Jaco_hole(j+1,j)= Jaco_hole(j+1,j)-0.5*dphidx-thermal/Deltax;
end
 
res_hole(1,1) = hole(1,1);
Jaco_hole(1,:) = 0.0;
Jaco_hole(1,1) = 1.0;    
res_hole(N,1) = elec(N,1) - Nd(N,1);
Jaco_hole(N,:) = 0.0;
Jaco_hole(N,N) = 1.0;
update_hole = Jaco_hole\(-res_hole);
hole= hole + update_hole;
figure(2);
plot (x, elec/1e6,'color', 'r');
hold on;
plot (x, hole/1e6,'color', 'b');
hold on;
