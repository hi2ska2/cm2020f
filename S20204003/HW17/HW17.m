gate_max=0;
voltage_step=0.05;
Max_VD=1.1;   % Targe drain voltage
i_VD=Max_VD/0.05+1;
Max_Vg=1.1;   % Target gate voltage
i_Vg=Max_Vg/0.05+1;
total_current=zeros(i_VD,i_Vg);


for gate_max=0:i_Vg-1
    
    clearvars -except total_current gate_max i_VD i_Vg
    
    % Parameters
    q=1.602192e-19;   % Elementary Charge, C
    e_si=11.7; e_ox=3.9;  e_mean=(e_si+e_ox)/2;  % Relative Permmitivity and mean Permitivity
    e0=8.8541878128e-12;   % Vacuum permittivity, F/m
    k_B= 1.380662e-23;  T=300;             % Boltzmann constant, J/K
    m0=9.109534e-31;
    nint=1.075e+16;  % Intrinsic carrier density, /m^3
    h=6.62607015e-34; hbar=h/(2*pi);           % Planck constant, J s
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
    dz= (bulk+2*tox)/(Nzz-1);
    
    % Doping 정보
    Nbody=0;   % Body Doping, 1/m^3, Polar : (+) for Donor , (-) for Acceptor
    Nd=1e+26;    % Drain, Source n-type Doping, 1/m^3
    
    
    %index 정보
    interface1=round(tox/dz)+1; interface2=round((tox+bulk)/dz)+1;   %interface index
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
    dy2=1 /dz^2;

    
    %행렬 구성
    Vbarrier=0.33374;
    Vg=0;
    VD=0;
    
    % For Jacobbian except for charges
    
    
    
    for gate_bias=0:gate_max
        
        Vg=(0.05*gate_bias);
        
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
        phi(interface1:interface2,Source)=thermal*asinh(Nd/(2*nint));
        phi(interface1:interface2,Drain)=thermal*asinh(Nd/(2*nint));
        phi=reshape(transpose(phi),[Nzz*Nxx,1]);
        
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
                            res(ix,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_si*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_ox*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/(2*e0);
                        end
                        
                    elseif iy==interface2  %Si-Oxide interface
                        
                        if ix==1+Nxx*(iy-1)   %neuman boundary
                            res(ix,1)=0;
                        elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                            res(ix,1)=0;
                        else
                            Jaco(ix,ix)=Jaco(ix,ix)-q*(nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/thermal/(2*e0);
                            res(ix,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_ox*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_si*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/(2*e0);
                        end
                        
                    elseif interface1<iy && iy<interface2   %Silicon
                        
                        if ix==1+Nxx*(iy-1) %neuman boundary
                            res(ix,1)=0;
                        elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                            res(ix,1)=0;
                        else
                            Jaco(ix,ix)=Jaco(ix,ix)-q*(nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/thermal/e0;
                            res(ix,1)=dx2*e_si*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_si*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-nint*exp(phi(ix,1)/thermal)+nint*exp(-phi(ix,1)/thermal))/e0;
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
            if max(abs(update_phi(:,Newton1)))<1e-15
                break;
            end
            
        end
        disp(sprintf('Gate voltage:%d  Drain voltage:%d \n', Vg , VD));
    end
    
    elec=zeros(Nxx*Nzz,1);
    elec((interface1-1)*Nxx+1:interface2*Nxx,1)=nint*exp(phi((interface1-1)*Nxx+1:interface2*Nxx,1)/thermal);
    
    hole=zeros(Nxx*Nzz,1);
    hole((interface1-1)*Nxx+1:interface2*Nxx,1)=nint*exp(-phi((interface1-1)*Nxx+1:interface2*Nxx,1)/thermal);
    
    
    A=sparse(3*Nxx*Nzz,3*Nxx*Nzz);
    
    
    for iy=1:Nzz    %y방향을 행 index로 하고 x방향을 열의 index로 하였습니다(좌표평면과 비슷)
        for ix=1+(iy-1)*Nxx:Nxx+(iy-1)*Nxx    %x축에 대한 행렬 구성이 끝나면 그 다음 y index에 대해 계산합니다
            if iy==1   %Oxide
                if ix>=Gi1 && ix<=Gi2
                    A(3*ix-2,3*ix-2)=1;
                elseif ix==1+Nxx*(iy-1)  %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox/2-dy2*e_ox/2; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox/2;
                elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox/2-dy2*e_ox/2; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox/2;
                else
                    A(3*ix-2,3*ix-2)=-dx2*e_ox-dy2*e_ox; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox;
                end
            elseif iy==Nzz   %Oxide
                if ix>=Gi1+Nxx*(iy-1) && ix<=Gi2+Nxx*(iy-1)
                    A(3*ix-2,3*ix-2)=1;
                elseif ix==1+Nxx*(iy-1)  %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox/2-dy2*e_ox/2; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox/2;
                elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox/2-dy2*e_ox/2; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox/2;
                else
                    A(3*ix-2,3*ix-2)=-dx2*e_ox-dy2*e_ox; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox/2; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox;
                end
            elseif iy==interface1     %Oxide-Si interface
                if ix==1+Nxx*(iy-1)   %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                else
                    A(3*ix-2,3*ix-2)=-2*dx2*e_mean-2*dy2*e_mean; A(3*ix-2,3*(ix-1)-2)=dx2*e_mean; A(3*ix-2,3*(ix+1)-2)=dx2*e_mean; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_si; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox;
                end
            elseif iy==interface2  %Si-Oxide interface
                if ix==1+Nxx*(iy-1)   %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                else
                    A(3*ix-2,3*ix-2)=-2*dx2*e_mean-2*dy2*e_mean; A(3*ix-2,3*(ix-1)-2)=dx2*e_mean; A(3*ix-2,3*(ix+1)-2)=dx2*e_mean; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_si;
                end
            elseif interface1<iy && iy<interface2   %Silicon
                if ix==1+Nxx*(iy-1) %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                elseif ix==Nxx+Nxx*(iy-1) %neuman boundary
                    A(3*ix-2,3*ix-2)=1;
                else
                    A(3*ix-2,3*ix-2)=-2*dx2*e_si-2*dy2*e_si; A(3*ix-2,3*(ix+1)-2)=dx2*e_si; A(3*ix-2,3*(ix-1)-2)=dx2*e_si; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_si; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_si;
                end
            elseif iy<interface1 || iy>interface2    %Oxide
                if ix==1+Nxx*(iy-1)  %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox-dy2*e_ox; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox/2; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox/2;
                elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                    A(3*ix-2,3*ix-2)=-dx2*e_ox-dy2*e_ox; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox/2; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox/2;
                else
                    A(3*ix-2,3*ix-2)=-2*dx2*e_ox-2*dy2*e_ox; A(3*ix-2,3*(ix+1)-2)=dx2*e_ox; A(3*ix-2,3*(ix-1)-2)=dx2*e_ox; A(3*ix-2,3*(ix+Nxx)-2)=dy2*e_ox; A(3*ix-2,3*(ix-Nxx)-2)=dy2*e_ox;
                end
            end
        end
    end
    for bias_drain=1:i_VD
        VD = (bias_drain-1)*0.05;
        
        for Newton2=1:100
            Jaco=sparse(3*Nxx*Nzz,1);
            res=zeros(3*Nxx*Nzz,1);
            Jaco=A;
            
            for iy=1:Nzz    %y방향을 행 index로 하고 x방향을 열의 index로 하였습니다(좌표평면과 비슷)
                for ix=1+(iy-1)*Nxx:Nxx+(iy-1)*Nxx    %x축에 대한 행렬 구성이 끝나면 그 다음 y index에 대해 계산합니다
                    
                    if iy==1
                        if ix>=Gi1 && ix<=Gi2      %Dirchlet Boundary
                            res(3*ix-2,1)=phi(ix,1)-Vbarrier-Vg;
                        elseif ix==1+Nxx*(iy-1)  %neuman boundary
                            res(3*ix-2,1)=0.5*dx2*e_ox*(phi(ix+1,1)-phi(ix,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                        elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                            res(3*ix-2,1)=-0.5*dx2*e_ox*(phi(ix,1)-phi(ix-1,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                        else
                            res(3*ix-2,1)=dx2*0.5*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_ox*(phi(ix+Nxx,1)-phi(ix,1));
                        end
                        
                    elseif iy==Nzz
                        if ix>=Gi1+Nxx*(iy-1) && ix<=Gi2+Nxx*(iy-1)  %Dirchlet Boundary
                            res(3*ix-2,1)=phi(ix,1)-Vbarrier-Vg;
                        elseif ix==1+Nxx*(iy-1)  %neuman boundary
                            res(3*ix-2,1)=0.5*dx2*e_ox*(phi(ix+1,1)-phi(ix,1))-0.5*dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                        elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                            res(3*ix-2,1)=-0.5*dx2*e_ox*(phi(ix,1)-phi(ix-1,1))-0.5*dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                        else
                            res(3*ix-2,1)=dx2*0.5*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))-dy2*e_ox*(phi(ix,1)-phi(ix-Nxx,1));
                        end
                        
                    elseif iy==interface1     %Oxide-Si interface
                        
                        if ix==1+Nxx*(iy-1)   %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint));
                        elseif ix==Nxx+Nxx*(iy-1) %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint))-VD;
                        else
                            Jaco(3*ix-2,3*ix-1)=-q/(e0*2);
                            Jaco(3*ix-2,3*ix)=q/(e0*2);
                            res(3*ix-2,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_si*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_ox*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-elec(ix,1)+hole(ix,1))/(2*e0);
                        end
                        
                    elseif iy==interface2  %Si-Oxide interface
                        
                        if ix==1+Nxx*(iy-1)   %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint));
                        elseif ix==Nxx+Nxx*(iy-1) %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint))-VD;
                        else
                            Jaco(3*ix-2,3*ix-1)=-q/(e0*2);
                            Jaco(3*ix-2,3*ix)=q/(e0*2);
                            res(3*ix-2,1)=dx2*e_mean*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*(e_ox*phi(ix+Nxx,1)-2*e_mean*phi(ix,1)+e_si*phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-elec(ix,1)+hole(ix,1))/(2*e0);
                        end
                        
                    elseif interface1<iy && iy<interface2   %Silicon
                        
                        if ix==1+Nxx*(iy-1) %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint));
                        elseif ix==Nxx+Nxx*(iy-1) %Dirichlet
                            res(3*ix-2,1)=phi(ix,1)-thermal*asinh(Nd/(2*nint))-VD;
                        else
                            Jaco(3*ix-2,3*ix-1)=-q/e0;
                            Jaco(3*ix-2,3*ix)=q/e0;
                            res(3*ix-2,1)=dx2*e_si*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_si*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1))+q*(Ndop(ix-Nxx*(iy-1),1)-elec(ix,1)+hole(ix,1))/e0;
                        end
                        
                    elseif iy<interface1 || iy>interface2    %Oxide
                        if ix==1+Nxx*(iy-1)  %neuman boundary
                            res(3*ix-2,1)=dx2*e_ox*(phi(ix+1,1)-phi(ix,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                        elseif ix==Nxx+Nxx*(iy-1)   %neuman boundary
                            res(3*ix-2,1)=-dx2*e_ox*(phi(ix,1)-phi(ix-1,1))+0.5*dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                        else
                            res(3*ix-2,1)=dx2*e_ox*(phi(ix+1,1)-2*phi(ix,1)+phi(ix-1,1))+dy2*e_ox*(phi(ix+Nxx,1)-2*phi(ix,1)+phi(ix-Nxx,1));
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for electron
                    
                    if iy>interface1-1 && iy<interface2+1
                        
                        if ix>1+(iy-1)*Nxx  && ix< Nxx+(iy-1)*Nxx
                            
                            x1=(phi(ix+1,1)-phi(ix,1))/thermal;
                            x2=(phi(ix,1)-phi(ix+1,1))/thermal;
                            x3=(phi(ix,1)-phi(ix-1,1))/thermal;
                            x4=(phi(ix-1,1)-phi(ix,1))/thermal;
                            
                            if abs((phi(ix+1,1)-phi(ix,1))/thermal)<0.02502 || abs((phi(ix,1)-phi(ix-1,1))/thermal)<0.02502
                                
                                Bern_P1 = ( 1.0-(x1)/2.0+(x1)^2/12.0*(1.0-(x1)^2/60.0*(1.0-(x1)^2/42.0)) ) ;
                                Bern_N1 = ( 1.0-(x2)/2.0+(x2)^2/12.0*(1.0-(x2)^2/60.0*(1.0-(x2)^2/42.0)) ) ;
                                Bern_P2 = ( 1.0-(x3)/2.0+(x3)^2/12.0*(1.0-(x3)^2/60.0*(1.0-(x3)^2/42.0)) ) ;
                                Bern_N2 = ( 1.0-(x4)/2.0+(x4)^2/12.0*(1.0-(x4)^2/60.0*(1.0-(x4)^2/42.0)) ) ;
                                
                                Deri_Bern_P1_phi1 = (-0.5 + (x1)/6.0*(1.0-(x1)^2/30.0*(1.0-(x1)^2/28.0)) )/thermal;
                                Deri_Bern_N1_phi1 = -(-0.5 + (x2)/6.0*(1.0-(x2)^2/30.0*(1.0-(x2)^2/28.0)) )/thermal;
                                Deri_Bern_P2_phi3 = -(-0.5 + (x3)/6.0*(1.0-(x3)^2/30.0*(1.0-(x3)^2/28.0)) )/thermal;
                                Deri_Bern_N2_phi3 = (-0.5 + (x4)/6.0*(1.0-(x4)^2/30.0*(1.0-(x4)^2/28.0)) )/thermal;
                                
                            elseif abs((phi(ix+1,1)-phi(ix,1))/thermal)<0.15 || abs((phi(ix,1)-phi(ix-1,1))/thermal)<0.15
                                
                                Bern_P1 = ( 1.0-(x1)/2.0+(x1)^2/12.0*(1.0-(x1)^2/60.0*(1.0-(x1)^2/42.0*(1-(x1)^2/40*(1-0.02525252525252525252525*(x1)^2)))));
                                Bern_N1 = ( 1.0-(x2)/2.0+(x2)^2/12.0*(1.0-(x2)^2/60.0*(1.0-(x2)^2/42.0*(1-(x2)^2/40*(1-0.02525252525252525252525*(x2)^2)))));
                                Bern_P2 = ( 1.0-(x3)/2.0+(x3)^2/12.0*(1.0-(x3)^2/60.0*(1.0-(x3)^2/42.0*(1-(x3)^2/40*(1-0.02525252525252525252525*(x3)^2)))));
                                Bern_N2 = ( 1.0-(x4)/2.0+(x4)^2/12.0*(1.0-(x4)^2/60.0*(1.0-(x4)^2/42.0*(1-(x4)^2/40*(1-0.02525252525252525252525*(x4)^2)))));
                                
                                Deri_Bern_P1_phi1 = (-0.5 + (x1)/6.0*(1.0-(x1)^2/30.0*(1.0-(x1)^2/28.0*(1-(x1)^2/30*(1-0.03156565656565656565657*(x1)^2)))))/thermal;
                                Deri_Bern_N1_phi1 = -(-0.5 + (x2)/6.0*(1.0-(x2)^2/30.0*(1.0-(x2)^2/28.0*(1-(x2)^2/30*(1-0.03156565656565656565657*(x2)^2)))))/thermal;
                                Deri_Bern_P2_phi3 = -(-0.5 + (x3)/6.0*(1.0-(x3)^2/30.0*(1.0-(x3)^2/28.0*(1-(x3)^2/30*(1-0.03156565656565656565657*(x3)^2)))))/thermal;
                                Deri_Bern_N2_phi3 = (-0.5 + (x4)/6.0*(1.0-(x4)^2/30.0*(1.0-(x4)^2/28.0*(1-(x4)^2/30*(1-0.03156565656565656565657*(x4)^2)))))/thermal;
                            elseif abs((phi(ix+1,1)-phi(ix,1))/thermal)>150.01 || abs((phi(ix,1)-phi(ix-1,1))/thermal)>150.01
                                
                                Bern_P1 = x1*exp(-x1);
                                Bern_N1 = x2*exp(-x2);
                                Bern_P2 = x3*exp(-x3);
                                Bern_N2 = x4*exp(-x4);
                                
                                Deri_Bern_P1_phi1 = (exp(-x1)-x1*exp(-x1))/thermal;
                                Deri_Bern_N1_phi1 = -(exp(-x2)-x2*exp(-x2))/thermal;
                                Deri_Bern_P2_phi3 = -(exp(-x3)-x3*exp(-x3))/thermal;
                                Deri_Bern_N2_phi3 = (exp(-x4)-x4*exp(-x4))/thermal;
                            else
                                
                                Bern_P1 = x1/(exp(x1)-1);
                                Bern_N1 = x2/(exp(x2)-1);
                                Bern_P2 = x3/(exp(x3)-1);
                                Bern_N2 = x4/(exp(x4)-1);
                                
                                Deri_Bern_P1_phi1=(1/(exp(x1)-1)-Bern_P1*(1/(exp(x1)-1)+1))/thermal;
                                Deri_Bern_N1_phi1=-(1/(exp(x2)-1)-Bern_N1*(1/(exp(x2)-1)+1))/thermal;
                                Deri_Bern_P2_phi3=-(1/(exp(x3)-1)-Bern_P2*(1/(exp(x3)-1)+1))/thermal;
                                Deri_Bern_N2_phi3=(1/(exp(x4)-1)-Bern_N2*(1/(exp(x4)-1)+1))/thermal;
                            end
                            
                            Deri_Bern_P1_phi2=-Deri_Bern_P1_phi1;
                            Deri_Bern_N1_phi2=-Deri_Bern_N1_phi1;
                            Deri_Bern_P2_phi2=-Deri_Bern_P2_phi3;
                            Deri_Bern_N2_phi2=-Deri_Bern_N2_phi3;
                            
                            res(3*ix-1,1) =  elec(ix+1,1)*(Bern_P1) - elec(ix,1)*(Bern_N1) - (elec(ix,1)*(Bern_P2) - elec(ix-1,1)*(Bern_N2));
                            
                            Jaco(3*ix-1,3*(ix+1)-1) = Bern_P1 ;
                            Jaco(3*ix-1,3*(ix+1)-2)= elec(ix+1,1)*Deri_Bern_P1_phi1 - elec(ix,1)*Deri_Bern_N1_phi1;
                            Jaco(3*ix-1,3*ix-1) =  - Bern_N1 - Bern_P2;
                            Jaco(3*ix-1,3*ix-2)= elec(ix+1,1)*Deri_Bern_P1_phi2 - elec(ix,1)*Deri_Bern_N1_phi2 - ( elec(ix,1)*Deri_Bern_P2_phi2 - elec(ix-1,1)*Deri_Bern_N2_phi2) ;
                            Jaco(3*ix-1,3*(ix-1)-1) = Bern_N2 ;
                            Jaco(3*ix-1,3*(ix-1)-2)= - (elec(ix,1)*Deri_Bern_P2_phi3 - elec(ix-1,1)*Deri_Bern_N2_phi3);
                            
                        else
                            res(3*ix-1,1) = elec(ix,1) - Nd;
                            Jaco(3*ix-1,:) = 0;
                            Jaco(3*ix-1,3*ix-1) = 1;
                        end
                        
                    else
                        Jaco(3*ix-1,3*ix-1)=1;
                        
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for holes
                    if iy>interface1-1 && iy<interface2+1
                        
                        if ix>1+(iy-1)*Nxx  && ix< Nxx+(iy-1)*Nxx
                            
                            x1=(phi(ix+1,1)-phi(ix,1))/thermal;
                            x2=(phi(ix,1)-phi(ix+1,1))/thermal;
                            x3=(phi(ix,1)-phi(ix-1,1))/thermal;
                            x4=(phi(ix-1,1)-phi(ix,1))/thermal;
                            
                            if abs((phi(ix+1,1)-phi(ix,1))/thermal)<0.02502 || abs((phi(ix,1)-phi(ix-1,1))/thermal)<0.02502
                                
                                Bern_P1 = ( 1.0-(x1)/2.0+(x1)^2/12.0*(1.0-(x1)^2/60.0*(1.0-(x1)^2/42.0)) ) ;
                                Bern_N1 = ( 1.0-(x2)/2.0+(x2)^2/12.0*(1.0-(x2)^2/60.0*(1.0-(x2)^2/42.0)) ) ;
                                Bern_P2 = ( 1.0-(x3)/2.0+(x3)^2/12.0*(1.0-(x3)^2/60.0*(1.0-(x3)^2/42.0)) ) ;
                                Bern_N2 = ( 1.0-(x4)/2.0+(x4)^2/12.0*(1.0-(x4)^2/60.0*(1.0-(x4)^2/42.0)) ) ;
                                
                                Deri_Bern_P1_phi1 = (-0.5 + (x1)/6.0*(1.0-(x1)^2/30.0*(1.0-(x1)^2/28.0)) )/thermal;
                                Deri_Bern_N1_phi1 = -(-0.5 + (x2)/6.0*(1.0-(x2)^2/30.0*(1.0-(x2)^2/28.0)) )/thermal;
                                Deri_Bern_P2_phi3 = -(-0.5 + (x3)/6.0*(1.0-(x3)^2/30.0*(1.0-(x3)^2/28.0)) )/thermal;
                                Deri_Bern_N2_phi3 = (-0.5 + (x4)/6.0*(1.0-(x4)^2/30.0*(1.0-(x4)^2/28.0)) )/thermal;
                                
                            elseif abs((phi(ix+1,1)-phi(ix,1))/thermal)<0.15 || abs((phi(ix,1)-phi(ix-1,1))/thermal)<0.15
                                
                                Bern_P1 = ( 1.0-(x1)/2.0+(x1)^2/12.0*(1.0-(x1)^2/60.0*(1.0-(x1)^2/42.0*(1-(x1)^2/40*(1-0.02525252525252525252525*(x1)^2)))));
                                Bern_N1 = ( 1.0-(x2)/2.0+(x2)^2/12.0*(1.0-(x2)^2/60.0*(1.0-(x2)^2/42.0*(1-(x2)^2/40*(1-0.02525252525252525252525*(x2)^2)))));
                                Bern_P2 = ( 1.0-(x3)/2.0+(x3)^2/12.0*(1.0-(x3)^2/60.0*(1.0-(x3)^2/42.0*(1-(x3)^2/40*(1-0.02525252525252525252525*(x3)^2)))));
                                Bern_N2 = ( 1.0-(x4)/2.0+(x4)^2/12.0*(1.0-(x4)^2/60.0*(1.0-(x4)^2/42.0*(1-(x4)^2/40*(1-0.02525252525252525252525*(x4)^2)))));
                                
                                Deri_Bern_P1_phi1 = (-0.5 + (x1)/6.0*(1.0-(x1)^2/30.0*(1.0-(x1)^2/28.0*(1-(x1)^2/30*(1-0.03156565656565656565657*(x1)^2)))))/thermal;
                                Deri_Bern_N1_phi1 = -(-0.5 + (x2)/6.0*(1.0-(x2)^2/30.0*(1.0-(x2)^2/28.0*(1-(x2)^2/30*(1-0.03156565656565656565657*(x2)^2)))))/thermal;
                                Deri_Bern_P2_phi3 = -(-0.5 + (x3)/6.0*(1.0-(x3)^2/30.0*(1.0-(x3)^2/28.0*(1-(x3)^2/30*(1-0.03156565656565656565657*(x3)^2)))))/thermal;
                                Deri_Bern_N2_phi3 = (-0.5 + (x4)/6.0*(1.0-(x4)^2/30.0*(1.0-(x4)^2/28.0*(1-(x4)^2/30*(1-0.03156565656565656565657*(x4)^2)))))/thermal;
                                
                            elseif abs((phi(ix+1,1)-phi(ix,1))/thermal)>150.01 || abs((phi(ix,1)-phi(ix-1,1))/thermal)>150.01
                                
                                Bern_P1 = x1*exp(-x1);
                                Bern_N1 = x2*exp(-x2);
                                Bern_P2 = x3*exp(-x3);
                                Bern_N2 = x4*exp(-x4);
                                
                                Deri_Bern_P1_phi1 = (exp(-x1)-x1*exp(-x1))/thermal;
                                Deri_Bern_N1_phi1 = -(exp(-x2)-x2*exp(-x2))/thermal;
                                Deri_Bern_P2_phi3 = -(exp(-x3)-x3*exp(-x3))/thermal;
                                Deri_Bern_N2_phi3 = (exp(-x4)-x4*exp(-x4))/thermal;
                                
                            else
                                
                                Bern_P1 = x1/(exp(x1)-1);
                                Bern_N1 = x2/(exp(x2)-1);
                                Bern_P2 = x3/(exp(x3)-1);
                                Bern_N2 = x4/(exp(x4)-1);
                                
                                Deri_Bern_P1_phi1=(1/(exp(x1)-1)-Bern_P1*(1/(exp(x1)-1)+1))/thermal;
                                Deri_Bern_N1_phi1=-(1/(exp(x2)-1)-Bern_N1*(1/(exp(x2)-1)+1))/thermal;
                                Deri_Bern_P2_phi3=-(1/(exp(x3)-1)-Bern_P2*(1/(exp(x3)-1)+1))/thermal;
                                Deri_Bern_N2_phi3=(1/(exp(x4)-1)-Bern_N2*(1/(exp(x4)-1)+1))/thermal;
                                
                            end
                            
                            Deri_Bern_P1_phi2=-Deri_Bern_P1_phi1;
                            Deri_Bern_N1_phi2=-Deri_Bern_N1_phi1;
                            Deri_Bern_P2_phi2=-Deri_Bern_P2_phi3;
                            Deri_Bern_N2_phi2=-Deri_Bern_N2_phi3;
                            
                            res(3*ix,1) =  hole(ix+1,1)*(Bern_N1) - hole(ix,1)*(Bern_P1) - (hole(ix,1)*(Bern_N2) - hole(ix-1,1)*(Bern_P2));
                            
                            Jaco(3*ix,3*(ix+1)) = Bern_N1 ;
                            Jaco(3*ix,3*(ix+1)-2)= hole(ix+1,1)*Deri_Bern_N1_phi1 - hole(ix,1)*Deri_Bern_P1_phi1;
                            Jaco(3*ix,3*ix) =  - Bern_P1 - Bern_N2;
                            Jaco(3*ix,3*ix-2)= hole(ix+1,1)*Deri_Bern_N1_phi2 - hole(ix,1)*Deri_Bern_P1_phi2 - ( hole(ix,1)*Deri_Bern_N2_phi2 - hole(ix-1,1)*Deri_Bern_P2_phi2) ;
                            Jaco(3*ix,3*(ix-1)) = Bern_P2 ;
                            Jaco(3*ix,3*(ix-1)-2)= - (hole(ix,1)*Deri_Bern_N2_phi3 - hole(ix-1,1)*Deri_Bern_P2_phi3);
                            
                        else
                            res(3*ix,1) = hole(ix,1) - (nint^2/Nd);
                            Jaco(3*ix,:) = 0;
                            Jaco(3*ix,3*ix) = 1;
                        end
                        
                    else
                        Jaco(3*ix,3*ix)=1;
                        
                    end
                    
                end
                
            end
            
            
            Cvector= zeros(3*Nxx*Nzz,1);
            Cvector(1:3:3*Nxx*Nzz-2,1) = thermal;
            Cvector(2:3:3*Nxx*Nzz-1,1)=Nd;
            Cvector(3:3:3*Nxx*Nzz,1)=nint^2/Nd;
            Cmatrix = spdiags(Cvector,0,3*Nxx*Nzz,3*Nxx*Nzz);
            Jaco_scaled = Jaco * Cmatrix;
            Rvector = 1./sum(abs(Jaco_scaled),2);
            Rmatrix = spdiags(Rvector,0,3*Nxx*Nzz,3*Nxx*Nzz);
            Jaco_scaled = Rmatrix* Jaco_scaled;
            res_scaled = Rmatrix *res;
            
            
            update_scaled=Jaco_scaled \ (-res_scaled);
            update_vector= Cmatrix* update_scaled;

            
            phi(:,1)=phi(:,1)+update_vector(1:3:3*Nxx*Nzz-2,1);
            elec(:,1)=elec(:,1)+update_vector(2:3:3*Nxx*Nzz-1,1);
            hole(:,1)=hole(:,1)+update_vector(3:3:3*Nxx*Nzz,1);
            
            if max(abs(update_vector(1:3:3*Nxx*Nzz-2,1)))<1e-15
                break;
            end
            
        end
        %    save_phi(:,:,iVg)=transpose(reshape(phi,[Nxx,Nzz]));  % 게이트 전압에 따른 Potential save
        disp(sprintf('Gate voltage:%d  Drain voltage:%d \n', Vg , VD));
        
        phi_c=transpose(reshape(phi,[Nxx,Nzz]));
        elec_c=transpose(reshape(elec,[Nxx,Nzz]));
        hole_c=transpose(reshape(hole,[Nxx,Nzz]));
        
        Current=zeros(Nzz,1);
        for iz=interface1:interface2
            if iz ~=interface1 && iz ~=interface2
                Current(iz,1)=q*1417*( elec_c(iz,Nxx) * (phi_c(iz,Nxx)-phi_c(iz,Nxx-1))/dx-thermal*((elec_c(iz,Nxx)-elec_c(iz,Nxx-1))/dx) )/1e+4*dz;
            elseif iz==interface1 || interface2
                Current(iz,1)=q*1417*( elec_c(iz,Nxx) * (phi_c(iz,Nxx)-phi_c(iz,Nxx-1))/dx-thermal*((elec_c(iz,Nxx)-elec_c(iz,Nxx-1))/dx) )/1e+4/2*dz;
            end
            total_current(bias_drain,gate_max+1)=sum(Current);
        end
        
        
    end
    
    
    phi=transpose(reshape(phi,[Nxx,Nzz]));
    elec=transpose(reshape(elec,[Nxx,Nzz]));
    hole=transpose(reshape(hole,[Nxx,Nzz]));
    
    
    plot(0.05*(1:bias_drain),total_current);
end

