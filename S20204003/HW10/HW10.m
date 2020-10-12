% Parameters
q=1.602192e-19;   % Elementary Charge, C
e_si=11.7; e_ox=3.9;  e_mean=(e_si+e_ox)/2;  % Relative Permmitivity and mean Permitivity
e0=8.854e-12;   % Vacuum permittivity, F/m
k_B= 1.380662e-23;  T=300;             % Boltzmann constant, J/K
m0=9.109534e-31;
nint=1.075e+16;  % Intrinsic carrier density, /m^3
h=6.626e-34; hbar=h/(2*pi);           % Planck constant, J s
% Device 길이 정보
GLg=10e-9;      %Gate length, m
CLg=GLg;
DLg=10e-9;      %Drain length, m
SLg=DLg;      %Source length, m
FLg=CLg+DLg+SLg;      %Channel length, m
tox=1e-9;             % Oxide thickness, m
bulk=5e-9;                   %Silicon thickness, m

% Mesh 정보
Nxx=61;          % X axis mesh
Nzz=36;          % Y axis mesh
dx= FLg/(Nxx-1);
dy= (bulk+2*tox)/(Nzz-1);


% Doping 정보
Nbody=0;   % Body Doping, 1/m^3, Polar : (+) for Donor , (-) for Acceptor
Nd=1e+26;    % Drain, Source n-type Doping, 1/m^3


%index 정보
interface1=round(tox/dy)+1; interface2=round((tox+bulk)/dy)+1;   %interface index
Gi1=(Nxx+1)/2-GLg/dx/2; Gi2=(Nxx+1)/2+GLg/dx/2;        %Gate index
Source=int16(SLg/dx+1); Drain=int16((SLg+CLg)/dx+1); Ndop=zeros(Nxx,1);  %Source and Drain index
Ndop(:,1)=Nbody;                                        %Channel Doping
Ndop(1:Source-1,1)=Nd;                                  %Source Doping
Ndop(Drain+1:Nxx,1)=Nd;                                 %Drain Doping
Ndop(Source,1)=(Nd+Nbody)/2;                               %Source/Channel interface doping charge
Ndop(Drain,1)=(Nd+Nbody)/2;                                %Channel/Drain interface doping charge

thermal=k_B*T/q;
A=sparse(Nxx*Nzz,Nxx*Nzz);
Jaco=sparse(Nxx*Nzz,Nxx*Nzz);
res=zeros(Nxx*Nzz,1);
dx2=1 /dx^2;
dy2=1 /dy^2;
f=zeros(Nxx*Nzz,1);
N1=interface2-interface1-1;

%행렬 구성
Vbarrier=0.33;
Elec=zeros(1,Nxx);
eDensity=zeros(N1,Nxx);
% For Jacobbian except for charges

for iy=1:Nzz    %y방향을 행 index로 하고 x방향을 열의 index로 하였습니다(좌표평면과 비슷)
    for ix=1+(iy-1)*Nxx:Nxx+(iy-1)*Nxx    %x축에 대한 행렬 구성이 끝나면 그 다음 y index에 대해 계산합니다
        if iy==1   %Oxide
            if ix>=Gi1 && ix<=Gi2
                A(ix,ix)=1;
            elseif ix==1+Nxx*(iy-1)  %neuman boundary
                A(ix,ix)=-dx2*e_ox/2-dy2*e_ox/2; A(ix,ix+1)=dx2*e_ox/2; A(ix,ix+Nxx)=dy2*e_ox/2;
            elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                A(ix,ix)=-dx2*e_ox/2-dy2*e_ox/2; A(ix,ix-1)=dx2*e_ox/2; A(ix,ix+Nxx)=dy2*e_ox/2;
            else
                A(ix,ix)=-dx2*e_ox-dy2*e_ox; A(ix,ix+1)=dx2*e_ox/2; A(ix,ix-1)=dx2*e_ox/2; A(ix,ix+Nxx)=dy2*e_ox;
            end
        elseif iy==Nzz   %Oxide
            if ix>=Gi1+Nxx*(iy-1) && ix<=Gi2+Nxx*(iy-1)
                A(ix,ix)=1;
            elseif ix==1+Nxx*(iy-1)  %neuman boundary
                A(ix,ix)=-dx2*e_ox/2-dy2*e_ox/2; A(ix,ix+1)=dx2*e_ox/2; A(ix,ix-Nxx)=dy2*e_ox/2;
            elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                A(ix,ix)=-dx2*e_ox/2-dy2*e_ox/2; A(ix,ix-1)=dx2*e_ox/2; A(ix,ix-Nxx)=dy2*e_ox/2;
            else
                A(ix,ix)=-dx2*e_ox-dy2*e_ox; A(ix,ix+1)=dx2*e_ox/2; A(ix,ix-1)=dx2*e_ox/2; A(ix,ix-Nxx)=dy2*e_ox;
            end
        elseif iy==interface1     %Oxide-Si interface
            if ix==1+Nxx*(iy-1)   %neuman boundary
                A(ix,ix)=1;
            elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                A(ix,ix)=1;
            else
                A(ix,ix)=-2*dx2*e_mean-2*dy2*e_mean; A(ix,ix-1)=dx2*e_mean; A(ix,ix+1)=dx2*e_mean; A(ix,ix+Nxx)=dy2*e_si; A(ix,ix-Nxx)=dy2*e_ox;
            end
        elseif iy==interface2  %Si-Oxide interface
            if ix==1+Nxx*(iy-1)   %neuman boundary
                A(ix,ix)=1;
            elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                A(ix,ix)=1;
            else
                A(ix,ix)=-2*dx2*e_mean-2*dy2*e_mean; A(ix,ix-1)=dx2*e_mean; A(ix,ix+1)=dx2*e_mean; A(ix,ix+Nxx)=dy2*e_ox; A(ix,ix-Nxx)=dy2*e_si;
            end
        elseif interface1<iy && iy<interface2   %Silicon
            if ix==1+Nxx*(iy-1) %neuman boundary
                A(ix,ix)=1;
            elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                A(ix,ix)=1;
            else
                A(ix,ix)=-2*dx2*e_si-2*dy2*e_si; A(ix,ix+1)=dx2*e_si; A(ix,ix-1)=dx2*e_si; A(ix,ix+Nxx)=dy2*e_si; A(ix,ix-Nxx)=dy2*e_si;
            end
        elseif iy<interface1 || iy>interface2    %Oxide
            if ix==1+Nxx*(iy-1)  %neuman boundary
                A(ix,ix)=-dx2*e_ox-dy2*e_ox; A(ix,ix+1)=dx2*e_ox; A(ix,ix+Nxx)=dy2*e_ox/2; A(ix,ix-Nxx)=dy2*e_ox/2;
            elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                A(ix,ix)=-dx2*e_ox-dy2*e_ox; A(ix,ix-1)=dx2*e_ox; A(ix,ix+Nxx)=dy2*e_ox/2; A(ix,ix-Nxx)=dy2*e_ox/2;
            else
                A(ix,ix)=-2*dx2*e_ox-2*dy2*e_ox; A(ix,ix+1)=dx2*e_ox; A(ix,ix-1)=dx2*e_ox; A(ix,ix+Nxx)=dy2*e_ox; A(ix,ix-Nxx)=dy2*e_ox;
            end
        end
    end
end

% Initial Guess
phi=zeros(Nzz,Nxx);
phi(interface1:interface2,1:Source-1)=thermal*asinh(Nd/(2*nint));
phi(interface1:interface2,Drain+1:Nxx)=thermal*asinh(Nd/(2*nint));
phi(interface1:interface2,Source)=thermal*asinh(Nd/(4*nint));
phi(interface1:interface2,Drain)=thermal*asinh(Nd/(4*nint));
phi=reshape(transpose(phi),[Nzz*Nxx,1]);

for iVg=1:12
    Vg=(iVg-1)/10;
    for Newton1=1:500
        Jaco=A;
        for iy=1:Nzz    %y방향을 행 index로 하고 x방향을 열의 index로 하였습니다(좌표평면과 비슷)
            
            for ix=1+(iy-1)*Nxx:Nxx+(iy-1)*Nxx    %x축에 대한 행렬 구성이 끝나면 그 다음 y index에 대해 계산합니다
                
                if iy==1
                    if ix>=Gi1 && ix<=Gi2      %Dirchlet Boundary
                        res(ix,1)=phi(ix,1)-Vbarrier-Vg;
                    elseif ix==1+Nxx*(iy-1)  %neuman boundary
                        res(ix,1)=0.5*dx2*e_ox*(phi(ix+1,1)-phi(ix,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                    elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                        res(ix,1)=-0.5*dx2*e_ox*(phi(ix,1)-phi(ix-1,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                    else
                        res(ix,1)=dx2*0.5*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                    end
                    
                elseif iy==Nzz
                    if ix>=Gi1+Nxx*(iy-1) && ix<=Gi2+Nxx*(iy-1)  %Dirchlet Boundary
                        res(ix,1)=phi(ix,1)-Vbarrier-Vg;
                    elseif ix==1+Nxx*(iy-1)  %neuman boundary
                        res(ix,1)=0.5*dx2*e_ox*(phi(ix+1,1)-phi(ix,1))-0.5*dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                    elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                        res(ix,1)=-0.5*dx2*e_ox*(phi(ix,1)-phi(ix-1,1))-0.5*dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                    else
                        res(ix,1)=dx2*0.5*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))-dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                    end
                    
                elseif iy==interface1     %Oxide-Si interface
                    
                    if ix==1+Nxx*(iy-1)   %neuman boundary
                        res(ix,1)=0;
                    elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                        res(ix,1)=0;
                    else
                        Jaco(ix,ix)=Jaco(ix,ix)-q*(nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/thermal/(2*e0);
                        res(ix,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_si*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_ox*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)-nint*exp(-phi(ix,1)/thermal))/(2*e0);
                    end
                    
                elseif iy==interface2  %Si-Oxide interface
                    
                    if ix==1+Nxx*(iy-1)   %neuman boundary
                        res(ix,1)=0;
                    elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                        res(ix,1)=0;
                    else
                        Jaco(ix,ix)=Jaco(ix,ix)-q*(nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/thermal/(2*e0);
                        res(ix,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_ox*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_si*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)-nint*exp(-phi(ix,1)/thermal))/(2*e0);
                    end
                    
                elseif interface1<iy && iy<interface2   %Silicon
                    
                    if ix==1+Nxx*(iy-1) %neuman boundary
                        res(ix,1)=0;
                    elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                        res(ix,1)=0;
                    else
                        Jaco(ix,ix)=Jaco(ix,ix)-q*(nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/thermal/e0;
                        res(ix,1)=dx2*e_si*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_si*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)-nint*exp(-phi(ix,1)/thermal))/e0;
                    end
                    
                elseif iy<interface1 || iy>interface2    %Oxide
                    if ix==1+Nxx*(iy-1)  %neuman boundary
                        res(ix,1)=dx2*e_ox*(phi(ix+1,1)-phi(ix,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                    elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                        res(ix,1)=-dx2*e_ox*(phi(ix,1)-phi(ix-1,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                    else
                        res(ix,1)=dx2*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                    end  
                end
            end
        end
        
        update_phi(:,Newton1) = (Jaco) \ (-res);
        phi = phi + update_phi(:,Newton1);
        
        error(Newton1,iVg)=max(abs(update_phi(:,Newton1)));
        if max(abs(update_phi(:,Newton1)))<1e-15
            break;
        end
        
    end
    save_phi(:,:,iVg)=transpose(reshape(phi,[Nxx,Nzz]));  % 게이트 전압에 따른 Potential save
end

phi=transpose(reshape(phi,[Nxx,Nzz]));


% Plot
x=0:dx*1e9:30;
y=0:dy*1e9:7;
surf(x,y,phi);
xlabel('Transport Direction (nm)');
ylabel('Confined Direction (nm)');
zlabel('Electrostatic Potential (V)');
