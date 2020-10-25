clear;
%%%%
q = 1.602e-19;
eps0 = 8.8542e-12;
k_B = 1.380662e-23;
T = 300;
V_T = k_B*T/q;
h = 6.626e-34;
hbar = h/(2*pi);
m0 = 9.1e-31;
Lx = 100e-9; Ly = 100e-9; Lz =5e-9;
mxx = 0.19; myy = 0.19; mzz = 0.91;
nmax = 10;
coef_2 = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Ec_Ei = 0.561004;


Si = 5e-9;
Nz = 121;
Ox = 0.5e-9;
delta_z = (Si+2*Ox)/(Nz-1);
eps_si = 11.7; eps_ox = 3.9;
Nacc = 1e24;
ni = 1.075e16;
coef = delta_z^2*q/eps0;
Vg = 0:0.02:1;
z = -Ox-Si/2:delta_z:Ox+Si/2;
phi = zeros(Nz,length(Vg));
ED = zeros(Nz,length(Vg));



interface1 = round((Nz)*Ox/(Si+2*Ox));
interface2 = round((Nz)*(Ox+Si)/(Si+2*Ox));
iteration = 100;
iter = 1:1:iteration;
error = zeros(iteration,length(Vg));
for k = 1:length(Vg)
    if (k>1)
        phi(:,k) = phi(:,k-1);
    end
    res = zeros(Nz,1);
    Jaco = zeros(Nz,Nz);
    phi(1,k) = 0.33374+Vg(1,k);
    phi(Nz,k) = 0.33374+Vg(1,k);
    res(1,1) = phi(1,k) - 0.33374-Vg(1,k);
    res(Nz,1) = phi(Nz,k) - 0.33374-Vg(1,k);
    Jaco(1,1) = 1; Jaco(Nz,Nz) = 1;
    for kk = 1:iteration
        if (kk==1)
            for j=1:20 %% Poisson iteration
                for i=2:interface1
                    if (i==interface1)
                        res(i,1) = eps_si*phi(i+1,k)-(eps_si+eps_ox)*phi(i,k)+eps_ox*phi(i-1,k)-coef*(Nacc+ni*exp(phi(i,k)/V_T))*0.5;
                        Jaco(i,i+1) = eps_si; Jaco(i,i) = -(eps_si+eps_ox)-coef*ni*exp(phi(i,k)/V_T)/V_T*0.5; Jaco(i,i-1) = eps_ox;
                    else
                        res(i,1) = eps_ox*phi(i+1,k)-2*eps_ox*phi(i,k)+eps_ox*phi(i-1,k);
                        Jaco(i,i-1) = eps_ox; Jaco(i,i) = -2*eps_ox; Jaco(i,i+1) = eps_ox;
                    end
                    
                end
                for i=interface2:Nz-1
                    if (i==interface2)
                        res(i,1) = eps_si*phi(i-1,k)-(eps_si+eps_ox)*phi(i,k)+eps_ox*phi(i+1,k)-coef*(Nacc+ni*exp(phi(i,k)/V_T))*0.5;
                        Jaco(i,i-1) = eps_si; Jaco(i,i) = -(eps_si+eps_ox)-coef*ni*exp(phi(i,k)/V_T)/V_T*0.5; Jaco(i,i+1) = eps_ox;
                    else
                        res(i,1) = eps_ox*phi(i+1,k)-2*eps_ox*phi(i,k)+eps_ox*phi(i-1,k);
                        Jaco(i,i-1) = eps_ox; Jaco(i,i) = -2*eps_ox; Jaco(i,i+1) = eps_ox;
                    end
                end
                for i=interface1+1:interface2-1
                    res(i,1) = eps_si*phi(i+1,k)-2*eps_si*phi(i,k)+eps_si*phi(i-1,k)-coef*(Nacc+ni*exp(phi(i,k)/V_T));
                    Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si-coef*ni*exp(phi(i,k)/V_T)/V_T; Jaco(i,i+1) = eps_si;
                end
                update = inv(Jaco)*(-res);
                phi(:,k) = phi(:,k) +update;
            end           
        else
            for j=1:1 %% Poisson iteration
%                 phi(:,k) = V_T*log(ED(:,k)/ni);
                for i=2:interface1
                    if (i==interface1)
                        res(i,1) = eps_si*phi(i+1,k)-(eps_si+eps_ox)*phi(i,k)+eps_ox*phi(i-1,k)-coef*(Nacc+ED(i,k))*0.5;
                        Jaco(i,i+1) = eps_si; Jaco(i,i) = -(eps_si+eps_ox)-coef*ni*exp(phi(i,k)/V_T)/V_T*0.5; Jaco(i,i-1) = eps_ox;
                    else
                        res(i,1) = eps_ox*phi(i+1,k)-2*eps_ox*phi(i,k)+eps_ox*phi(i-1,k);
                        Jaco(i,i-1) = eps_ox; Jaco(i,i) = -2*eps_ox; Jaco(i,i+1) = eps_ox;
                    end
                    
                end
                for i=interface2:Nz-1
                    if (i==interface2)
                        res(i,1) = eps_si*phi(i-1,k)-(eps_si+eps_ox)*phi(i,k)+eps_ox*phi(i+1,k)-coef*(Nacc+ED(i,k))*0.5;
                        Jaco(i,i-1) = eps_si; Jaco(i,i) = -(eps_si+eps_ox)-coef*ni*exp(phi(i,k)/V_T)/V_T*0.5; Jaco(i,i+1) = eps_ox;
                    else
                        res(i,1) = eps_ox*phi(i+1,k)-2*eps_ox*phi(i,k)+eps_ox*phi(i-1,k);
                        Jaco(i,i-1) = eps_ox; Jaco(i,i) = -2*eps_ox; Jaco(i,i+1) = eps_ox;
                    end
                end
                for i=interface1+1:interface2-1
                    res(i,1) = eps_si*phi(i+1,k)-2*eps_si*phi(i,k)+eps_si*phi(i-1,k)-coef*(ED(i,k)+Nacc);
                    Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si-coef*ni*exp(phi(i,k)/V_T)/V_T; Jaco(i,i+1) = eps_si;
                end
                update = inv(Jaco)*(-res);
                phi(:,k) = phi(:,k) +update;
            end  
        end
        error(kk,k) = max(update);
        %%%        
        V = q*Ec_Ei-q*phi(:,k);
        Nbulk = interface2-interface1-1;
        Hamil = zeros(Nbulk,Nbulk);
        for ii=2:Nbulk-1
            Hamil(ii,ii-1) = 1;
            Hamil(ii,ii) = -2;
            Hamil(ii,ii+1) = 1;
        end
        Hamil(1,1) = -2; Hamil(1,2) = 1; Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
        
        for ii = 1:Nbulk
            Hamil(ii,ii) = Hamil(ii,ii) - 2*mzz*m0*(delta_z/hbar)^2*V(ii+interface1,1);
            %         Hamil(ii,ii) = Hamil(ii,ii) -2*mass(3)*m0*(Deltaz/hbar)^2*V(ii+interface1,i);
        end
        
        [eigenvector,eigenvalues] = eig(Hamil);
        scaledEz = diag(eigenvalues)/(-2*mzz*m0*(delta_z/hbar)^2);
        [sortedEz,sortedIndex] = sort(scaledEz);
        
        elec = zeros(Nz,1);
        TotalNumber = 0;
        for n=1:nmax
            Ez = sortedEz(n,1);
            wavefunction2 = eigenvector(:,sortedIndex(n)).^2;
            normalization = sum(wavefunction2)*delta_z;
            wavefunction2 = wavefunction2/normalization;
            SubbandNumber = coef_2*log(1+exp(-(Ez)/(k_B*T)));
            TotalNumber = TotalNumber+SubbandNumber;
            elec(interface1+1:interface2-1,1) = elec(interface1+1:interface2-1,1) + 1/(Lx*Ly)*wavefunction2*SubbandNumber;
        end
        ED(:,k) = elec;
    end
end

figure(1)
plot(z,ED*1e-6);
xlabel('position z (nm)');
ylabel('Charge density (/cm^3)');
figure(2)
plot(iter,error);
xlabel('Iteration');
ylabel('Error (Max)');
figure(3)
semilogy(iter,error);
xlabel('Iteration');
ylabel('Error (Max)');








