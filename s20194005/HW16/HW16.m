clear;
q = 1.602192e-19; 
eps0 = 8.854187e-12; 
Deltax = 0.2e-9;
m0 = 9.109383e-31 ;
N = 301; 
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
interface1 = 51;
interface2 = 251;
Nd= zeros(N,1);
Nd(1:interface1,1) = 5e23;
Nd(interface1+1:interface2-1,1) = 2e21;
Nd(interface2:N,1) = 5e23;
res = zeros (N,1);
Jaco = sparse (N,N);
phi = zeros (N,1);
phi(:,1)=thermal*log(Nd(:,1)/ni);
res2 = zeros (2*N,1);
Jaco2 = sparse(2*N,2*N);   
mobility = 2000e-4;
Area = 1e-12; 
Current=zeros(11,1);
for newton=1:20

    res(1,1)= phi(1,1) -thermal*log(Nd(1,1)/ni);
    Jaco(1,1) =1;
    for j=2:N-1
        res(j,1)= eps_si*(phi(j+1,1)-2*phi(j,1)+phi(j-1,1))-coef*(-Nd(j,1)+ni*exp(phi(j,1)/thermal));
        Jaco(j,j+1) = eps_si;
        Jaco(j,j) = -2*eps_si-coef*ni*exp(phi(j,1)/thermal)/thermal;
        Jaco(j,j-1) = eps_si;
    end
    res(N,1) = phi(N,1) - thermal*log(Nd(1,1)/ni);
    Jaco(N,N) = 1;
    update = Jaco\(-res);
    phi = phi + update;
end 
figure(1);
elec1 = zeros(N,1);
elec1 = ni*exp(phi/thermal);
plot(x,elec1/1e6,'color', 'r');
hold on

%continuity equation
elec2 = zeros(N,1);
elec2 = ni*exp(phi/thermal); 
for bias= 1:11
    Voltage=(bias-1)*0.05;
for nt = 1:20
    res2(1,1) = phi(1,1) - thermal*log(Nd(1,1)/ni);
    Jaco2(1,1) = 1;
     for j = 2:N-1
        res2(2*j-1,1) = eps_si*(phi(j+1,1)-2*phi(j,1)+phi(j-1,1))+coef*(Nd(j,1)-elec2(j,1));
        Jaco2(2*j-1,2*j+1) = eps_si;
        Jaco2(2*j-1,2*j-1) = -2*eps_si;
        Jaco2(2*j-1,2*j-3) = eps_si;
        Jaco2(2*j-1,2*j) = -coef;
     end
    res2(2*N-1,1) = phi(N,1) - thermal*log(Nd(N,1)/ni) - Voltage;
    Jaco2(2*N-1,2*N-1) = 1;
for j=1:N-1
    n_av=0.5*(elec2(j+1,1)+elec2(j,1));
    dphidx= (phi(j+1,1)-phi(j,1))/Deltax;
    delecdx=(elec2(j+1,1)-elec2(j,1))/Deltax;
    Jn= n_av * dphidx- thermal *delecdx;
    res2(2*j,1)=res2(2*j,1)+ Jn;
    Jaco2(2*j,2*j+2)= Jaco2(2*j,2*j+2)+0.5*dphidx-thermal/Deltax;
    Jaco2(2*j,2*j)= Jaco2(2*j,2*j)+0.5*dphidx+thermal/Deltax;
    Jaco2(2*j,2*j+1)= Jaco2(2*j,2*j+1)+n_av/Deltax;
    Jaco2(2*j,2*j-1)= Jaco2(2*j,2*j-1)-n_av/Deltax;
    res2(2*j+2,1)= res2(2*j+2,1)-Jn;
    Jaco2(2*j+2,2*j+2)=Jaco2(2*j+2,2*j+2)-0.5*dphidx+thermal/Deltax;
    Jaco2(2*j+2,2*j)= Jaco2(2*j+2,2*j)-0.5*dphidx-thermal/Deltax;
    Jaco2(2*j+2,2*j+1)= Jaco2(2*j+2,2*j+1)-n_av/Deltax;
    Jaco2(2*j+2,2*j-1)= Jaco2(2*j+2,2*j-1)+n_av/Deltax;
end

 %Scaling
res2(2,1) = elec2(1,1)-Nd(1,1);
Jaco2(2,:) = 0.0;
Jaco2(2,2) = 1.0;    
res2(2*N,1) = elec2(N,1) - Nd(N,1);
Jaco2(2*N,:) = 0.0;
Jaco2(2*N,2*N) = 1.0;
Cvector= zeros(2*N,1);
Cvector(1:2:2*N-1,1) = thermal;
Cvector(2:2:2*N,1)=max(abs(Nd));
Cmatrix = spdiags(Cvector,0,2*N,2*N);
Jaco_scaled = Jaco2 * Cmatrix;
Rvector = 1./sum(abs(Jaco_scaled),2);
Rmatrix = spdiags(Rvector,0,2*N,2*N);
Jaco_scaled = Rmatrix* Jaco_scaled;
res_scaled = Rmatrix *res2;
update_scaled=Jaco_scaled \ (-res_scaled);
update= Cmatrix* update_scaled;
phi = phi + update(1:2:2*N-1,1);
elec2 = elec2 + update(2:2:2*N,1);
end
Current(bias,1)=q*mobility*(elec2(N,1)*(phi(N,1)-phi(N-1,1))/Deltax-thermal*(elec2(N,1)-elec2(N-1,1))/Deltax)/1e+4*Area;
end
figure(2);
plot (x, elec2/1e6,'color', 'b');
hold on;
figure(3);
plot (0:0.05:0.5,Current/1e6);
hold on;
