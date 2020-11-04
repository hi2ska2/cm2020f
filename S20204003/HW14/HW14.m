q=1.602192e-19;     % Elementary charge, C
Ptype=300e-9;  % P type region length, m
Ntype=300e-9; % N type region length, m
N=601;   % Discretization
dx=(Ptype+Ntype)/(N-1);
nint=1.075e16;   % intrinsic carrier density, m-3
Ndop=1e23;   % N type doping concentration, m-3
Pdop=1e23;   % P type doping concentration, m-3
e1=11.7; e0=8.854187817e-12;   % Permittivity, F/m
k_B=1.380662e-23;  % Boltzmann constant, J/K
T=300;
jacob=zeros(N,N);
res=zeros(N,1);
interface=301;
thermal =k_B*T/q;
coeff=dx*dx*q/e0;


phi(1:interface,1) = thermal*asinh(-Pdop/2/nint);
phi(interface+1:N,1) = thermal*asinh(Ndop/2/nint);


for Newton=1:50
    for i=2:N-1
        res(1,1) = phi(1,1) - thermal*asinh(-Pdop/2/nint);
        jacob(1,1) = 1;
        if i<=interface
            res(i,1)=e1*(phi(i-1,1)-2*phi(i,1)+phi(i+1,1))+coeff*(-Pdop-nint*exp(phi(i,1)/thermal)+nint*exp(-phi(i,1)/thermal));
            jacob(i,i+1)=e1; jacob(i,i)=-2*e1-coeff*nint*(exp(phi(i,1)/thermal)+exp(-phi(i,1)/thermal))/thermal;
            jacob(i,i-1)=e1;
        elseif  i>interface
            res(i,1)=e1*(phi(i-1,1)-2*phi(i,1)+phi(i+1,1))+coeff*(Ndop-nint*exp(phi(i,1)/thermal)+nint*exp(-phi(i,1)/thermal));
            jacob(i,i+1)=e1; jacob(i,i)=-2*e1-coeff*nint*(exp(phi(i,1)/thermal)+exp(-phi(i,1)/thermal))/thermal;
            jacob(i,i-1)=e1;
        end
        res(N,1) = phi(N,1) - thermal*asinh(Ndop/2/nint);
        jacob(N,N) = 1;
        
    end
    
    delphi=jacob\(-res);
    savedelphi(:,Newton)=delphi;
    phi=phi+delphi;
end

%Plot
x=0:dx*1e+9:600;
elec = nint*exp(phi/thermal);
hole = nint*exp(-phi/thermal);

% Answer 1
elec_Poisson=elec/1e+6;
hole_Poisson=hole/1e+6;

plot(x,elec_Poisson,'r','linewidth',1.5)
hold on
plot(x,hole_Poisson,'k','linewidth',1.5)
hold on


%For electron continuity
res_elec = zeros(N,1);
jacob_elec = zeros(N,N);
jacob_elec(1,:) = 0;
jacob_elec(1,1) = 1;
res_elec(1,1) = elec(1,1) - nint^2/Pdop;
for ii=2:N-1
    n_av1 = 0.5*(elec(ii+1,1)+elec(ii,1));
    n_av2 = 0.5*(elec(ii,1)+elec(ii-1,1));
    dphidx1 = (phi(ii+1,1)-phi(ii,1))/dx;
    dphidx2 = (phi(ii,1)-phi(ii-1,1))/dx;
    delecdx1 = (elec(ii+1,1)-elec(ii,1))/dx;
    delecdx2 = (elec(ii,1)-elec(ii-1,1))/dx;
    res_elec(ii,1) =  n_av1 * dphidx1 - thermal * delecdx1 - n_av2 * dphidx2 + thermal * delecdx2 ;
    jacob_elec(ii,ii+1) = 0.5* dphidx1 -thermal/dx;
    jacob_elec(ii,ii) = 0.5* dphidx1 +thermal*2/dx -0.5*dphidx2;
    jacob_elec(ii,ii-1) = -0.5*dphidx2 - thermal/dx;
end
res_elec(N,1) = elec(N,1) - Ndop;
jacob_elec(N,:) = 0;
jacob_elec(N,N) = 1;

update_elec = jacob_elec \ (-res_elec);
elec = elec + update_elec;


%Hole continuity
res_hole = zeros(N,1);
jacob_hole = zeros(N,N);
jacob_hole(1,:) = 0;
jacob_hole(1,1) = 1;
res_hole(1,1) = hole(1,1) - Pdop;
for ii=2:N-1
    h_av1 = 0.5*(hole(ii+1,1)+hole(ii,1));
    h_av2 = 0.5*(hole(ii,1)+hole(ii-1,1));
    dphidx1 = (phi(ii+1,1)-phi(ii,1))/dx;
    dphidx2 = (phi(ii,1)-phi(ii-1,1))/dx;
    dholedx1 = (hole(ii+1,1)-hole(ii,1))/dx;
    dholedx2 = (hole(ii,1)-hole(ii-1,1))/dx;
    res_hole(ii,1) =  h_av1 * dphidx1 + thermal * dholedx1 - ( h_av2 * dphidx2 + thermal * dholedx2 );
    jacob_hole(ii,ii+1) = 0.5* dphidx1 + thermal/dx;
    jacob_hole(ii,ii) = 0.5* dphidx1 - thermal*2/dx -0.5*dphidx2;
    jacob_hole(ii,ii-1) = -0.5*dphidx2 + thermal/dx;
end
res_hole(N,1) = hole(N,1) - nint^2/Ndop;
jacob_hole(N,:) = 0;
jacob_hole(N,N) = 1;

update_hole = jacob_hole \ (-res_hole);
hole = hole + update_hole;


% Answer 2
elec_DD=elec/1e+6;
hole_DD=hole/1e+6;

%Plot
plot(x,elec_DD,'o','linewidth',1.1,'Displayname','Electron');
hold on;
plot(x,hole_DD,'or','linewidth',1.1,'Displayname','Hole');

xlabel('position(nm)')
ylabel('Electron Density(1/cm^3)')

