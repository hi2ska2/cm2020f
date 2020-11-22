q=1.602e-19;    % Elementary charge, C
e_si=11.7;
e0=8.854187817e-12;   % Vacuum permittivity, F/m
k_B=1.380662e-23;    % Boltzmann constant, J/K
T=300;                % Temperature, K
thermal=k_B*T/q; % Thermal voltage, V
HiDop=1e-8; %m  for N+ doped region.
LowDop=4e-8; %m  for N doped region.
N=301;
dx=(2*HiDop+LowDop)/(N-1);
nint=1.075e16;    % intrinsic carrier density, 1/m^3
Ndop1=5e+25;  %  Donor doped for N+ region, 1/m^3
Ndop2=2e+23;  % Donor doped for N region, 1/m^3
jacob=sparse(N,N);
res=zeros(N,1);
interface1=round(HiDop/dx+1);
interface2=round((HiDop+LowDop)/dx+1);
coeff=dx*dx*q/e0;

phi(1:interface1,1) = k_B*T/q*log(Ndop1/nint);
phi(interface1+1:interface2-1,1) = k_B*T/q*log(Ndop2/nint);
phi(interface2:N,1) = k_B*T/q*log(Ndop1/nint);

Mobility = 1417; % cm^2/V s
Area = 1e-4; % Area : 1 cm^2

% non-linear Poisson solver
for Newton=1:30
    for ii=2:N-1
        res(1,1) = phi(1,1) - thermal *log(Ndop1/nint);
        jacob(1,1) = 1;
        if ii<=interface1
            res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop1-nint*exp(phi(ii,1)/thermal));
            jacob(ii,ii+1)=e_si; jacob(ii,ii)=-2*e_si-coeff*nint*exp(phi(ii,1)/thermal)/thermal;
            jacob(ii,ii-1)=e_si;
        elseif  ii>interface1  && ii<interface2
            res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop2-nint*exp(phi(ii,1)/thermal));
            jacob(ii,ii+1)=e_si; jacob(ii,ii)=-2*e_si-coeff*nint*exp(phi(ii,1)/thermal)/thermal;
            jacob(ii,ii-1)=e_si;
        elseif  ii>=interface2
            res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop1-nint*exp(phi(ii,1)/thermal));
            jacob(ii,ii+1)=e_si; jacob(ii,ii)=-2*e_si-coeff*nint*exp(phi(ii,1)/thermal)/thermal;
            jacob(ii,ii-1)=e_si;
        end
        res(N,1) = phi(N,1) - k_B*T/q*log(Ndop1/nint);
        jacob(N,N) = 1;
        
    end
    
    
    delphi=jacob\(-res);
    phi=phi+delphi;
    
end

elec = nint*exp(phi/thermal);

res_elec = zeros(N,1);
jacob_elec = sparse(N,N);


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
    res_elec(1,1) = elec(1,1) - Ndop1;
    jacob_elec(1,:) = 0;
    jacob_elec(1,1) = 1;
    res_elec(N,1) = elec(N,1) - Ndop1;
    jacob_elec(N,:) = 0;
    jacob_elec(N,N) = 1;
end

update_elec = jacob_elec \ (-res_elec);
save(:,1) = update_elec;
elec = elec + update_elec;



Current=zeros(0.5/0.05+1,1);

for bias=0:10
    
    Voltage=(0.05*bias);
    
    for Newton=1:50
        Jaco=sparse(2*N,2*N);
        res=zeros(2*N,1);
        res(1,1) = phi(1,1) - k_B*T/q*log(Ndop1/nint);
        Jaco(1,1) = 1;
        for ii=2:N-1
            if ii<=interface1
                res(2*ii-1,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop1-elec(ii,1));
                Jaco(2*ii-1,2*ii+1)=e_si; Jaco(2*ii-1,2*ii-1)=-2*e_si; Jaco(2*ii-1,2*ii)=-coeff; Jaco(2*ii-1,2*(ii-1)-1)=e_si;
            elseif  ii>interface1  && ii<interface2
                res(2*ii-1,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop2-elec(ii,1));
                Jaco(2*ii-1,2*ii+1)=e_si; Jaco(2*ii-1,2*ii-1)=-2*e_si; Jaco(2*ii-1,2*ii)=-coeff; Jaco(2*ii-1,2*(ii-1)-1)=e_si;
            elseif  ii>=interface2
                res(2*ii-1,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(Ndop1-elec(ii,1));
                Jaco(2*ii-1,2*ii+1)=e_si; Jaco(2*ii-1,2*ii-1)=-2*e_si; Jaco(2*ii-1,2*ii)=-coeff; Jaco(2*ii-1,2*(ii-1)-1)=e_si;
            end
            res(2*N-1,1) = phi(N,1) - k_B*T/q*log(Ndop1/nint)- Voltage;
            Jaco(2*N-1,2*N-1) = 1;
            n_av1 = 0.5*(elec(ii+1,1)+elec(ii,1));
            n_av2 = 0.5*(elec(ii,1)+elec(ii-1,1));
            dphidx1 = (phi(ii+1,1)-phi(ii,1))/dx;
            dphidx2 = (phi(ii,1)-phi(ii-1,1))/dx;
            delecdx1 = (elec(ii+1,1)-elec(ii,1))/dx;
            delecdx2 = (elec(ii,1)-elec(ii-1,1))/dx;
            res(2*ii,1) =  n_av1 * dphidx1 - thermal * delecdx1 - n_av2 * dphidx2 + thermal * delecdx2 ;
            
            Jaco(2*ii,2*(ii+1)) = 0.5* dphidx1 -thermal/dx;
            Jaco(2*ii,2*ii+1)=n_av1/dx;
            Jaco(2*ii,2*ii) = 0.5* dphidx1 +thermal*2/dx -0.5*dphidx2;
            Jaco(2*ii,2*ii-1)=-n_av1/dx-n_av2/dx;
            Jaco(2*ii,2*(ii-1)) = -0.5*dphidx2 - thermal/dx;
            Jaco(2*ii,2*ii-3)=n_av2/dx;
        end
        res(2,1) = elec(1,1) - Ndop1;
        Jaco(2,:) = 0;
        Jaco(2,2) = 1;
        res(2*N,1) = elec(N,1) - Ndop1;
        Jaco(2*N,:) = 0;
        Jaco(2*N,2*N) = 1;
        
        Cvector= zeros(2*N,1);
        Cvector(1:2:2*N-1,1) = thermal;
        Cvector(2:2:2*N,1)=Ndop1;
        Cmatrix = spdiags(Cvector,0,2*N,2*N);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,2*N,2*N);
        Jaco_scaled = Rmatrix* Jaco_scaled;
        res_scaled = Rmatrix *res;
        
        update_scaled=Jaco_scaled \ (-res_scaled);
        update_vector= Cmatrix* update_scaled;
        
        phi(:,1)=phi(:,1)+update_vector(1:2:2*N-1,1);
        elec(:,1)=elec(:,1)+update_vector(2:2:2*N,1);
        
        save_update(:,Newton)=update_vector;
        maximum_update_phi(Newton,1)=max(abs(save_update(1:2:2*N-1,Newton)));
        maximum_update_elec(Newton,1)=max(abs(save_update(2:2:2*N,Newton)))/1e6;  % cm-3
        
    end
    % Current
    Current(bias+1,1)=q*Mobility*((elec(N,1)+elec(N-1,1))/2*(phi(N,1)-phi(N-1,1))/dx-thermal*(elec(N,1)-elec(N-1))/dx)/1e+4 * Area; %A/cm^2 * cm^2
    disp(sprintf('Current Voltage:%d \n', Voltage));
end

