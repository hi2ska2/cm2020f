clear;
q = 1.602192e-19; 
eps0 = 8.854187e-12; 
Delta_z = 0.1e-9;
m0 = 9.109383e-31 ;
N = 67; 
T=300.0;
hbar = 1.054571e-34;
k = 1.380662e-23;
ni = 1.075e16;
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e24; % 1e18/cm^3, we can change it to 1e23,24,25 
Lx=100e-9;
Ly=100e-9;
interface1 = 9;
interface2 = 59;
Applied_Voltage = 0:0.1:1;
coef = Delta_z*Delta_z*q/eps0;
thermal = k*T/q;
phi=zeros(N,11);
phi(:,1) = 0.33374;
electron_S = zeros(N,11);
Ndensity=zeros(N,11);
res= zeros(N,1);
electron_P = zeros(N,11);
    for newton = 1:30
        res = zeros(N,1);
        Jaco = sparse(N,N);
        res(1,1) = phi(1,1) - 0.33374;
        res(N,1) = phi(N,1) - 0.33374;
        Jaco(1,1) = 1;
        Jaco(N,N) = 1;
     
        for i = 2:66
            if      (i< interface1 || i> interface2)
                res(i,1) = eps_ox*phi(i+1,1)-2*eps_ox*phi(i,1)+eps_ox*phi(i-1,1)-coef*(Nacc+ni*exp(phi(i,1)/thermal));
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -2*eps_ox-coef*ni*exp(phi(i,1)/thermal)/thermal; 
                Jaco (i,i+1) = eps_ox;
            elseif  (i == interface1)
                res(i,1) = eps_si*phi(i+1,1)-(eps_si+eps_ox)*phi(i,1)+eps_ox*phi(i-1,1)-coef*(Nacc+ni*exp(phi(i,1)/thermal))*0.5;
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -(eps_si+eps_ox)-coef*ni*exp(phi(i,1)/thermal)/thermal*0.5; 
                Jaco(i,i+1) = eps_si;
            elseif  (i == interface2)
                res(i,1) = eps_ox*phi(i+1,1)-(eps_ox+eps_si)*phi(i,1)+eps_si*phi(i-1,1)-coef*(Nacc+ni*exp(phi(i,1)/thermal))*0.5;
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -(eps_ox+eps_si)-coef*ni*exp(phi(i,1)/thermal)/thermal*0.5; 
                Jaco(i,i+1) = eps_ox;
            else
                res(i,1) = eps_si*phi(i+1,1)-2*eps_si*phi(i,1)+eps_si*phi(i-1,1)-coef*(Nacc+ni*exp(phi(i,1)/thermal));
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -2*eps_si-coef*ni*exp(phi(i,1)/thermal)/thermal; 
                Jaco(i,i+1) = eps_si;
            end
        end
        
        update = Jaco\(-res);
        phi(:,1) = phi(:,1) + update;

    
z=0.1*[0:(N-1)];
for j=1:11
     if j>1
        phi(:,j) = phi(:,j-1);
     end
    V(:,j) = - q*phi(:,j)+q*0.56 ;
    electron_S(:,j) = 0;
for valley=1:3
    if valley==1
        mxx=0.91*m0; myy=0.19*m0;  mzz=0.19*m0; 
    elseif valley==2
        mxx=0.19*m0; myy=0.91*m0;  mzz=0.19*m0;
    else
        mxx=0.19*m0; myy=0.19*m0;  mzz=0.91*m0;
    end

        
    H= zeros(49,49);
    H(1,1)=-2-2*mzz/(hbar)^2*Delta_z^2*V(interface1+1,1); H(1,2)=1;
    H(49,49)=-2; H(49,48)=1;
    for i=2:48
        H(i,i-1)=1;
        H(i,i)=-2-2*mzz/(hbar)^2*Delta_z^2*V(interface1+1,j);
        H(i,i+1)=1;
    end
    
    [Eigenvectors,Eigenvalues] = eig(H);
    [S_Ez,S_Index]= sort(diag(Eigenvalues)/(-2*mzz*Delta_z^2)*hbar^2);
    Eigenvector_S=Eigenvectors(:,S_Index);
    normalize=zeros(49,49);
    for n=1:49
        Wavefunc=Eigenvector_S(:,n).^2;
        Normal=sum(Wavefunc*Delta_z);
        Wavefunc=Wavefunc/Normal;
        subbandNumber = 2*2*Lx*Ly/(2*pi)*sqrt(mxx*myy)/(hbar^2)*(k*T)*log(1+exp(-S_Ez(n,1)/(k*T)));
        electron_S(interface1+1:interface2-1,j) = electron_S(interface1+1:interface2-1,j) + 1/(Lx*Ly)*subbandNumber*Wavefunc;
    end
    end
        
   
    for newtons = 1:30
        res = zeros(N,1);
        Jaco = sparse(N,N);
        res(1,1) = phi(1,j) - 0.33374-0.1*(j-1);
        res(N,1) = phi(N,j) - 0.33374-0.1*(j-1);
        Jaco(1,1) = 1;
        Jaco(N,N) = 1;
        for i = 2:66
            if      (i< interface1 || i> interface2)
                res(i,1) = eps_ox*phi(i+1,j)-2*eps_ox*phi(i,j)+eps_ox*phi(i-1,j)-coef*(Nacc+electron_S(i,j));
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -2*eps_ox-coef*electron_S(i,j)/thermal; 
                Jaco(i,i+1) = eps_ox;
            elseif  (i == interface1)
                res(i,1) = eps_si*phi(i+1,j)-(eps_si+eps_ox)*phi(i,j)+eps_ox*phi(i-1,j)-coef*Nacc*0.5;
                Jaco(i,i-1) = eps_ox; 
                Jaco(i,i) = -(eps_si+eps_ox); 
                Jaco(i,i+1) = eps_si;
            elseif  (i == interface2)
                res(i,1) = eps_ox*phi(i+1,j)-(eps_ox+eps_si)*phi(i,j)+eps_si*phi(i-1,j)-coef*Nacc*0.5;
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -(eps_ox+eps_si); 
                Jaco(i,i+1) = eps_ox;
            else
                res(i,1) = eps_si*phi(i+1,j)-2*eps_si*phi(i,j)+eps_si*phi(i-1,j)-coef*(Nacc+electron_S(i,j));
                Jaco(i,i-1) = eps_si; 
                Jaco(i,i) = -2*eps_si-coef*electron_S(i,j)/thermal; 
                Jaco(i,i+1) = eps_si;
            end
        end

        update = Jaco\(-res);
        phi(:,j) = phi(:,j) + update;
        Max_P_Update = max(abs(update));
    
for m = interface1+1:interface2-1
        electron_S(m,j) =  electron_S(m,j)*exp(update(m,1)/(thermal));
end

    end
end
end

electron_S = electron_S/1e6;  
plot(z,electron_S);      

xlabel('position (nm)')
ylabel('Electron density (cm^-2)')
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1.0V');
