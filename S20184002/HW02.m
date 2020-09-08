clear;
%******
h = 6.63e-34/(2*pi); % reduced Planck's constant
m_e = 0.19*9.1e-31; % effective mass 
a = 5e-9; % width
N = 500; % discretization number N
delta_x = a/(N-1);
%******
for i=1:N-2
    v1(1,i) = -2;
end

for i=1:N-3
    v2(1,i) = 1;
end

A1 = diag(v1);
A2 = diag(v2,1);
A3 = diag(v2,-1);
A= A1+A2+A3;

[V,D] = eig(A);
M_1=D(N-2,N-2);
M_2=D(N-3,N-3);
M_3=D(N-4,N-4);

WaveFunction=(1/sqrt(delta_x)).*V;
Probability=WaveFunction.*WaveFunction;

E_eV_1 = M_1*h^2*(-1/delta_x^2)/(2*m_e)/(1.6e-19); % n=1 (numerical)
E_eV_2 = M_2*h^2*(-1/delta_x^2)/(2*m_e)/(1.6e-19); % n=2 (numerical)
E_eV_3 = M_3*h^2*(-1/delta_x^2)/(2*m_e)/(1.6e-19); % n=3 (numerical)

Analytic_E_eV_1 = h^2*pi^2/(2*m_e*a^2)/(1.6e-19); % n=1 (analytical)
Analytic_E_eV_2 = 4*h^2*pi^2/(2*m_e*a^2)/(1.6e-19); % n=2 (analytical)
Analytic_E_eV_3 = 9*h^2*pi^2/(2*m_e*a^2)/(1.6e-19); % n=3 (analytical)


