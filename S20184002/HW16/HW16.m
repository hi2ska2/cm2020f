clear;
q = 1.602192e-19;
eps0 = 8.8542e-12; % F/m
k_B = 1.380662e-23;
T = 300; % K
V_t = (k_B*T)/q;
X_Left = 0; X_Right = 600e-9;
XS = 100e-9;
XL = X_Right - X_Left;
Nx = 601;
deltaX = XL/(Nx-1);
x = X_Left:deltaX:X_Right;
x1 = (Nx-1)*(XS/XL)+1;
x2 = (Nx-1)*(XL-XS)/XL+1;
eps_si = 11.7; eps_ox = 3.9;
u_n = 0.135; % m^2/V*s
Area = 1e-12; % m^2
Npp = 5e23;
Np = 2e21;
Ndon = Npp*ones(Nx,1);
Ndon(x1+1:x2,1) = Np; 
ni = 1.075e16; % m^-3
coef = deltaX^2*q/eps0;
phi = zeros(Nx,1);
phi(1:Nx,1) = V_t*log(Ndon(1:Nx,1)/ni);

res = zeros(Nx,1);
Jaco = sparse(Nx,Nx);
Newton = 10;



for j=1:Newton
    res(1,1) = phi(1,1) - V_t*log(Ndon(1,1)/ni);
    res(Nx,1) = phi(Nx,1) - V_t*log(Ndon(Nx,1)/ni);
    Jaco(1,1) = 1; Jaco(Nx,Nx) = 1;   
    for i=2:Nx-1
        res(i,1) = eps_si*(phi(i+1,1)-2*phi(i,1)+phi(i-1,1))-coef*(-Ndon(i,1) + ni*exp(phi(i,1)/V_t) - ni*exp(-phi(i,1)/V_t));
        Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si - coef*ni*(exp(phi(i,1)/V_t) + exp(-phi(i,1)/V_t))/V_t; Jaco(i,i+1) = eps_si;
    end
    update = Jaco\(-res);
    phi = phi + update;
%     norm(update,inf);
end

elec = zeros(Nx,1);
elec = ni*exp(phi/V_t);

% plot (x,phi);
% figure(2)
% plot(x,elec);
elec_N = elec;
phi_N = phi;

BiasV = 0:0.05:0.5;
Current = zeros(length(BiasV),1);

for Bias = 1:length(BiasV)
    for newton=1:10
        res = zeros(2*Nx,1);
        Jaco = sparse(2*Nx,2*Nx);
        res(1,1) = phi(1,1) - V_t*log(Ndon(1,1)/ni);
        Jaco(1,1) = 1;
        
        for i=2:Nx-1
            res(2*i-1,1) = eps_si*(phi(i+1,1) - 2*phi(i,1) + phi(i-1,1)) + coef*(Ndon(i,1)-elec(i,1));
            Jaco(2*i-1,2*i+1) = eps_si;
            Jaco(2*i-1,2*i-1) = -2*eps_si;
            Jaco(2*i-1,2*i-3) = eps_si;
            Jaco(2*i-1,2*i) = -coef;
        end
        
        res(2*Nx-1,1) = phi(Nx,1) - V_t*log(Ndon(Nx,1)/ni) - BiasV(1,Bias);
        Jaco(2*Nx-1,2*Nx-1) = 1.0;
        
        for i=1:Nx-1
            n_av = 0.5*(elec(i+1,1)+elec(i,1));
            dphidx = (phi(i+1,1) - phi(i,1))/deltaX;
            delecdx = (elec(i+1,1) - elec(i,1))/deltaX;
            Jn = n_av * dphidx - V_t * delecdx;
            res(2*i,1) = res(2*i,1) + Jn;
            Jaco(2*i,2*i+2) = Jaco(2*i,2*i+2) + 0.5*dphidx - V_t / deltaX;
            Jaco(2*i,2*i) = Jaco(2*i,2*i) + 0.5*dphidx + V_t / deltaX;
            Jaco(2*i,2*i+1) = Jaco(2*i,2*i+1) + n_av / deltaX;
            Jaco(2*i,2*i-1) = Jaco(2*i,2*i-1) - n_av / deltaX;
            res(2*i+2,1) = res(2*i+2,1) - Jn;
            Jaco(2*i+2,2*i+2) = Jaco(2*i+2,2*i+2) - 0.5*dphidx + V_t /deltaX;
            Jaco(2*i+2,2*i) = Jaco(2*i+2,2*i) - 0.5*dphidx - V_t / deltaX;
            Jaco(2*i+2,2*i+1) = Jaco(2*i+2,2*i+1) - n_av / deltaX;
            Jaco(2*i+2,2*i-1) = Jaco(2*i+2,2*i-1) + n_av / deltaX;
        end
        
        res(2,1) = elec(1,1) - Ndon(1,1);
        Jaco(2,:) = 0;
        Jaco(2,2) = 1.0;
        res(2*Nx,1) = elec(Nx,1) - Ndon(Nx,1);
        Jaco(2*Nx,:) = 0;
        Jaco(2*Nx,2*Nx) = 1.0;
        
        Cvector = zeros(2*Nx,1);
        Cvector(1:2:2*Nx-1,1) = V_t;
        Cvector(2:2:2*Nx,1) = max(abs(Ndon));
        Cmatrix = spdiags(Cvector,0,2*Nx,2*Nx);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,2*Nx,2*Nx);
        Jaco_scaled = Rmatrix * Jaco_scaled;
        res_scaled = Rmatrix * res;
        update_scaled = Jaco_scaled \ (-res_scaled);
        update = Cmatrix * update_scaled;
        
        update = Jaco \ (-res);
        phi = phi + update(1:2:2*Nx-1,1);
        elec = elec + update(2:2:2*Nx,1);
        norm(update(1:2:2*Nx-1,1),inf);
    end
    Current(Bias,1) = Area*q*u_n*(0.5*(elec(Nx,1)+elec(Nx-1,1))*(phi(Nx,1)-phi(Nx-1,1))/deltaX - V_t*(elec(Nx,1)-elec(Nx-1,1))/deltaX);
end

% V = inv(Jaco_elec);
% figure(3)
% semilogy(x/1e-9,elec_N*1e-6,x/1e-9,elec*1e-6,'o');
% title('Electron');
% legend('Nonlinear','Continuity');
% xlabel('position x (nm)');
% ylabel('Carrier density (/cm^3)');
% 
% figure(4)
% semilogy(x/1e-9,abs(elec_N-elec)./elec);
% title('Error');
% 
% xlabel('position x (nm)');
% ylabel('Error');

plot(BiasV,1e3*Current,'-o');
title('Long');
xlabel('Bias Voltage (V)');
ylabel('I (mA)');

    





