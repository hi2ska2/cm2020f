clear;
q = 1.602192e-19;
eps0 = 8.8542e-12; % F/m
k_B = 1.380662e-23;
T = 300; % K
V_t = (k_B*T)/q;
X_Left = 0; X_Right = 600e-9;
XS = 300e-9;
XL = X_Right - X_Left;
Nx = 601;
deltaX = XL/(Nx-1);
x = X_Left:deltaX:X_Right;
x1 = (Nx-1)*(XS/XL)+1;
x2 = (Nx-1)*(XL-XS)/XL+1;
eps_si = 11.7; eps_ox = 3.9;
Ndon = 1e23*ones(Nx,1);
ni = 1.075e16; % m^-3
coef = deltaX^2*q/eps0;
phi = zeros(Nx,1);
phi(1:x1,1) = -V_t*log(Ndon(1:x1,1)/ni);
phi(x1+1:Nx,1) = V_t*log(Ndon(x1+1:Nx,1)/ni);
res = zeros(Nx,1);
Jaco = zeros(Nx,Nx);
Newton = 10;



for j=1:Newton
    res(1,1) = phi(1,1) + V_t*log(Ndon(1,1)/ni);
    res(Nx,1) = phi(Nx,1) - V_t*log(Ndon(Nx,1)/ni);
    Jaco(1,1) = 1; Jaco(Nx,Nx) = 1;   
    for i=2:x1
        res(i,1) = eps_si*(phi(i+1,1)-2*phi(i,1)+phi(i-1,1))-coef*(Ndon(i,1) + ni*exp(phi(i,1)/V_t) - ni*exp(-phi(i,1)/V_t));
        Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si - coef*ni*(exp(phi(i,1)/V_t) + exp(-phi(i,1)/V_t))/V_t; Jaco(i,i+1) = eps_si;
    end
    for i=x1+1:Nx-1
        res(i,1) = eps_si*(phi(i+1,1)-2*phi(i,1)+phi(i-1,1))-coef*(-Ndon(i,1) + ni*exp(phi(i,1)/V_t) - ni*exp(-phi(i,1)/V_t));
        Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si - coef*ni*(exp(phi(i,1)/V_t) + exp(-phi(i,1)/V_t))/V_t; Jaco(i,i+1) = eps_si;
    end
    update = Jaco\(-res);
    phi = phi + update;
%     norm(update,inf);
end

elec = zeros(Nx,1); hole = zeros(Nx,1);
elec = ni*exp(phi/V_t); hole = ni*exp(-phi/V_t);

plot (x,phi);
figure(2)
plot(x,elec,x,hole);

res_elec = zeros(Nx,1);
Jaco_elec = zeros(Nx,Nx);

res_elec(1,1) = elec(1,1) - ni^2/Ndon(1,1);
res_elec(Nx,1) = elec(Nx,1) - Ndon(Nx,1);
Jaco_elec(1,1) = 1; Jaco_elec(Nx,Nx) = 1;

for i=2:Nx-1
    n_av_b = 0.5*(elec(i,1) + elec(i-1,1));
    n_av_f = 0.5*(elec(i,1) + elec(i+1,1));
    dphidx_b = (phi(i,1) - phi(i-1,1))/deltaX;
    dphidx_f = (phi(i+1,1) - phi(i,1))/deltaX;
    delecdx_b = (elec(i,1) - elec(i-1,1))/deltaX;
    delecdx_f = (elec(i+1,1) - elec(i,1))/deltaX;
    res_elec(i,1) = (n_av_b*dphidx_f - V_t*delecdx_f) - (n_av_b*dphidx_b - V_t*delecdx_b);
    Jaco_elec(i,i-1) = -0.5*dphidx_b - V_t/deltaX;
    Jaco_elec(i,i) = 0.5*(dphidx_f - dphidx_b) - 2*V_t/deltaX; 
    Jaco_elec(i,i+1) = 0.5*dphidx_f - V_t/deltaX;
end
update_elec = Jaco_elec\(-res_elec);
elec_C = elec +update_elec;

%%% hole
res_hole = zeros(Nx,1);
Jaco_hole = zeros(Nx,Nx);

res_hole(1,1) = hole(1,1) - Ndon(1,1);
res_elec(Nx,1) = hole(Nx,1) - ni^2/Ndon(Nx,1);
Jaco_hole(1,1) = 1; Jaco_hole(Nx,Nx) = 1;

for i=2:Nx-1
    p_av_b = 0.5*(hole(i,1) + hole(i-1,1));
    p_av_f = 0.5*(hole(i,1) + hole(i+1,1));
    dphidx_b = (phi(i,1) - phi(i-1,1))/deltaX;
    dphidx_f = (phi(i+1,1) - phi(i,1))/deltaX;
    dholedx_b = (hole(i,1) - hole(i-1,1))/deltaX;
    dholedx_f = (hole(i+1,1) - hole(i,1))/deltaX;
    res_hole(i,1) = (p_av_b*dphidx_f + V_t*dholedx_f) - (p_av_b*dphidx_b + V_t*dholedx_b);
    Jaco_hole(i,i-1) = -0.5*dphidx_b + V_t/deltaX;
    Jaco_hole(i,i) = 0.5*(dphidx_f - dphidx_b) + 2*V_t/deltaX; 
    Jaco_hole(i,i+1) = 0.5*dphidx_f + V_t/deltaX;
end
update_hole = Jaco_hole\(-res_hole);
hole_C = hole + update_hole;



figure(3)
semilogy(x/1e-9,elec*1e-6,x/1e-9,elec_C*1e-6,'o');
title('Electron');
legend('Nonlinear','Continuity');
xlabel('position x (nm)');
ylabel('Carrier density (/cm^3)');
figure(4)
semilogy(x/1e-9,hole*1e-6,x/1e-9,hole_C*1e-6,'o');
title('Hole');
legend('Nonlinear','Continuity');
xlabel('position x (nm)');
ylabel('Carrier density (/cm^3)');
    





