clear;
%%%%%
V_T = 0.0259*1.6e-19; % Thermal voltage
C = zeros(1,18); % /m^3
for i=2:18
    C(1,i) = 10^(15+0.5*i);
end
ni = 1e16; % /m^3
N = 50; % Iteration number

% Analytic (+)
phi_A_p = V_T*asinh(C./(2*ni));
phi_A_eV_p = phi_A_p./1.6e-19;

% Analytic (-)
phi_A_m = V_T*asinh(-C./(2*ni));
phi_A_eV_m = phi_A_m./1.6e-19;

% Numerical (+)
phi_i_p = 0.35*1.6e-19;
J_p = C+ni*exp(-phi_i_p/V_T)-ni*exp(phi_i_p/V_T);
JJ_p = (ni*exp(-phi_i_p/V_T)+ni*exp(phi_i_p/V_T))*(-1/V_T);
phi_update_p = phi_i_p-J_p./JJ_p;

for i=1:N
    J_p = C+ni*exp(-phi_update_p/V_T)-ni*exp(phi_update_p/V_T);
    JJ_p = (ni*exp(-phi_update_p/V_T)+ni*exp(phi_update_p/V_T))*(-1/V_T);
    phi_update_p = phi_update_p-J_p./JJ_p;
end

% Numerical (-)
phi_i_m = -0.35*1.6e-19;
J_m = -C+ni*exp(-phi_i_m/V_T)-ni*exp(phi_i_m/V_T);
JJ_m = (ni*exp(-phi_i_m/V_T)+ni*exp(phi_i_m/V_T))*(-1/V_T);
phi_update_m = phi_i_m-J_m./JJ_m;

for i=1:N
    J_m = -C+ni*exp(-phi_update_m/V_T)-ni*exp(phi_update_m/V_T);
    JJ_m = (ni*exp(-phi_update_m/V_T)+ni*exp(phi_update_m/V_T))*(-1/V_T);
    phi_update_m = phi_update_m-J_m./JJ_m;
end

semilogx(C*1e-6,phi_A_eV_p,C*1e-6,phi_update_p./1.6e-19,'o',C*1e-6,phi_A_eV_m,C*1e-6,phi_update_m./1.6e-19,'s');
xlabel('Doping Concentration (cm^-3)') 
ylabel('Electrostatic Potential (eV)')
legend({'Positive (Analytical)','Positive (Numerical)','Negative (Analytical)','Negative (Analytical)'},'Location','southwest')



