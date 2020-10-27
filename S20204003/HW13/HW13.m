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


Electron=zeros(N,1);
Nz=zeros(N,1);

% initial guess
phi=zeros(N,1);
phi(:,1)=0.33374;


% 1. For Simple Poisson Solver


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
    Elec_save(:,index_Vg)=nint*exp(phi/thermal);
    Electron(interface1:interface2,1)=Elec_save(interface1:interface2,index_Vg); % Electron density, 1/cm^3
    eDensity_Poisson(:,index_Vg)=Electron/1e6;
end




% 2. For Self-Consistent Poisson-Schrodinger Solver

%intital guess
phi=zeros(N,1);
phi(:,1)=0.33374;

for Newton=1:40
    Jaco=sparse(N,N);
    res=zeros(N,1);
    for ii=1:N
        if ii==1 || ii==N
            res(ii,1)=phi(ii,1)-0.33374;
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


for index_Vg=1:11
    Vg=(index_Vg-1)*0.1;
    for Schrodinger=1:50
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
            for n=1:N1
                distribution=eigenvector_sorted(:,n).^2;
                Sum=sum(distribution*dz);
                normalize(:,n)=distribution/Sum;
            end
            
            for z=1:N1
                for n=1:nmax
                    subband(n,valley_type)=1*Lx*Ly/(2*pi)*md/(hbar^2)*k_B*T*log(1+exp((-Ezn(n,1))/(k_B*T)));
                    Elec_valley(z,valley_type)=Elec_valley(z,valley_type)+1/(Lx*Ly)*normalize(z,n)*subband(n,valley_type);
                end
            end
        end
        Nz(interface1+1:interface2-1,1)=2*2*sum(Elec_valley,2);
        
        for Newton1=1:30
            jacob=zeros(N,N);
            res=zeros(N,1);
            for ii=1:N
                if ii==1 || ii==N
                    res(ii,1)=phi(ii,1)-0.33374-Vg;
                    jacob(ii,ii)=1;
                elseif  interface1<ii && ii<interface2
                    res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(-Na-Nz(ii,1));
                    jacob(ii,ii+1)=e_si; jacob(ii,ii)=-2*e_si-coeff*Nz(ii,1)/thermal;
                    jacob(ii,ii-1)=e_si;
                elseif  ii<interface1 || ii>interface2
                    res(ii,1)=e_ox*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1));
                    jacob(ii,ii+1)=e_ox; jacob(ii,ii)=-2*e_ox;
                    jacob(ii,ii-1)=e_ox;
                elseif ii==interface1
                    res(ii,1)=e_si*(-phi(ii,1)+phi(ii+1,1))-e_ox*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-Nz(ii,1))/2;
                    jacob(ii,ii+1)=e_si;
                    jacob(ii,ii)=-e_si-e_ox-coeff*Nz(ii,1)/(2*thermal);
                    jacob(ii,ii-1)=e_ox;
                elseif ii==interface2
                    res(ii,1)=e_ox*(-phi(ii,1)+phi(ii+1,1))-e_si*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-Nz(ii,1))/2;
                    jacob(ii,ii+1)=e_ox;
                    jacob(ii,ii)=-e_si-e_ox-coeff*Nz(ii,1)/(2*thermal);
                    jacob(ii,ii-1)=e_si;
                end
            end
            update_phi2(:,Newton1)=jacob \ (-res);
            phi=phi+update_phi2(:,Newton1);
            Nz=Nz.*exp(update_phi2(:,Newton1)/thermal);
            
            if max(abs(update_phi2(:,Newton1)))<1e-15
                break;
            end
            
        end
        error(Schrodinger,index_Vg)=max(abs(update_phi2(:,1)));
        
        %   if max(abs(update_phi2(:,1)))<5e-12
        %       break;
        %   end
        
    end
    eDensity_Schrodinger(:,index_Vg)=Nz/1e6;
end
