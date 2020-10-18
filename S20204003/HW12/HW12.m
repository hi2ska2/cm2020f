% Parameters
q=1.602192e-19; % Elementary charge, C
m0=9.109e-31;  %electron mass, kg
h=6.626176e-34; hbar=h/(2*pi); %Planck constant, J s
k_B=1.38065e-23;     % Boltzmann constant, J/K
T=300;  % Temperature, K
e_si=11.7;  e_ox=3.9;  e0=8.854187817e-12; % Permittivity, F/m
nint=1.075e16; % Intrinsic carrier density, 1/m^3

% Device information
width=5e-9;                % Silicon thickness, m
tox=0.8e-9;           %Oxide thickness
N=67;
dz=(width+2*tox)/(N-1);
interface1=round(tox/dz)+1;  interface2=round((tox+width)/dz)+1;
N1=interface2-interface1-1;  %Silicon part discretization
Na=1e24;          % Doping concentration, 1/m^3
Lx=100e-9; Ly=100e-9;

nmax=N1;  % the number of subbands which are considered.
Jaco=sparse(N,N);
res=zeros(N,1);
coeff=dz*dz*q/e0;
thermal=k_B*T/q;


Electron=zeros(N,11);
Nz=zeros(N,11);

% initial guess
phi=zeros(N,1);
phi(:,1)=0.33374;


for index_Vg=1:11
    Vg=(index_Vg-1)*0.1; % Gate Voltage
    for Newton=1:40
        Jaco=sparse(N,N);
        res=zeros(N,1);
        for ii=1:N
            if ii==1 || ii==N
                res(ii,1)=phi(ii,1)-0.33374-Vg;
                Jaco(ii,ii)=1;
            elseif  interface1<ii && ii<interface2
                res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal));
                Jaco(ii,ii+1)=e_si; Jaco(ii,ii)=-2*e_si-coeff*nint*exp(phi(ii,1)/thermal)/thermal;
                Jaco(ii,ii-1)=e_si;
            elseif  ii<interface1 || ii>interface2
                res(ii,1)=e_ox*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1));
                Jaco(ii,ii+1)=e_ox; Jaco(ii,ii)=-2*e_ox; Jaco(ii,ii-1)=e_ox;
            elseif ii==interface1
                res(ii,1)=e_si*(-phi(ii,1)+phi(ii+1,1))-e_ox*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal))/2;
                Jaco(ii,ii+1)=e_si;
                Jaco(ii,ii)=-e_si-e_ox-coeff*nint*exp(phi(ii,1)/thermal)/(2*thermal);
                Jaco(ii,ii-1)=e_ox;
            elseif ii==interface2
                res(ii,1)=e_ox*(-phi(ii,1)+phi(ii+1,1))-e_si*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal))/2;
                Jaco(ii,ii+1)=e_ox;
                Jaco(ii,ii)=-e_si-e_ox-coeff*nint*exp(phi(ii,1)/thermal)/(2*thermal);
                Jaco(ii,ii-1)=e_si;
            end
        end
        delphi=Jaco\(-res);
        phi=phi+delphi;
        
        if max(abs(delphi))<1e-15    % Break
            break;
        end
    end
    Elec_save(:,index_Vg)=nint*exp(phi/thermal)/1e+6;
    Electron(interface1:interface2,index_Vg)=Elec_save(interface1:interface2,index_Vg); % Electron density, 1/cm^3
    
    Elec_valley=zeros(N1,3);
    ham=zeros(N1,N1);
    
    V=-q*phi+0.56*q;
    
    for valley_type=1:3
        
        if valley_type==1
            mzz=0.91*m0; mxx=0.19*m0; myy=0.19*m0; md=sqrt(mxx*myy);
        elseif valley_type==2
            mzz=0.19*m0; mxx=0.91*m0; myy=0.19*m0; md=sqrt(mxx*myy);
        elseif valley_type==3
            mzz=0.19*m0; mxx=0.19*m0; myy=0.91*m0; md=sqrt(mxx*myy);
        end
        
        for a=1:N1
            if a==1
                ham(a,a)=-2-2*mzz/(hbar)^2*dz*dz*V(interface1+a,1);
                ham(a,a+1)=1;
            elseif a==N1
                ham(a,a)=-2-2*mzz/(hbar)^2*dz*dz*V(interface1+a,1);
                ham(a,a-1)=1;
            else
                ham(a,a-1)=1;
                ham(a,a)=-2-2*mzz/(hbar)^2*dz*dz*V(interface1+a,1);
                ham(a,a+1)=1;
            end
        end
        [eigenvector,eigenvalue]=eig(ham);
        [Ezn,ind]=sort(diag(eigenvalue)/(-2*mzz*dz*dz)*hbar^2);
        eigenvector_sorted=eigenvector(:,ind);
        normalize=zeros(N1,N1);
        % Normalize eigenfunction
        for n=1:N1
            distribution=eigenvector_sorted(:,n).^2;
            Sum=sum(distribution*dz);
            normalize(:,n)=distribution/Sum;
        end
        
        for z=1:N1
            for n=1:nmax
                subband(n,valley_type)=Lx*Ly/(2*pi)*md/(hbar^2)*k_B*T*log(1+exp((-Ezn(n,1))/(k_B*T)));
                Elec_valley(z,valley_type)=Elec_valley(z,valley_type)+1/(Lx*Ly)*normalize(z,n)*subband(n,valley_type);
            end
        end
    end
    Nz(interface1+1:interface2-1,index_Vg)=2*2*sum(Elec_valley,2)/1e6;  %multiply 2 twice is for spin and valley degeneracies.
end

% Answer ( Row: at each position, Column: at each gate voltage)
eDensity_Poisson=Electron;
eDensity_Schrodinger=Nz;

