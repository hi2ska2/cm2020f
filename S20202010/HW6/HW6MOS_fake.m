%Assignment 6, MOS_fake
clear all;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, m
tsi = 5e-9; %silicon layer thickness, m

k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.075e16; % 1.075e10/cm^3 
T=300; %temp. 300K
Nacc = 1e24; %number density, m^-3, 10^18cm^-3
%prompt = 'What is the number N?';
C1 = q*Nacc*tsi/(2*eps0*eps_ox); %analytic constant A
C2 = q*Nacc/(2*eps0*eps_si); %analytic constant B
C3 = tsi/2+tox; %analytic constant b
C4 = q*Nacc*tox*tsi/(2*eps0*eps_ox)+q*Nacc*tsi*tsi*0.5*0.5/(2*eps0*eps_si); %analytic constant C
VG=zeros(11,1);
for kk=1:11
    VG(kk,1) = - (kk-1)/10.0;
end 

N = 67;

Deltax = (6.6e-9)/(N-1); %spacing
interface1 = 9; %at x =0.8nm
interface2 = 59; %at x =5.8nm
j1=interface1;
j2=interface2;
A=[];
A(1,1)=1;

y2=zeros(N,1);

for i=2:(N-1)
    if (i<j1) A(i,i-1)=eps_ox; A(i,i)=-2*eps_ox; A(i,i+1)=1*eps_ox;
    elseif (i==j1) A(i, i-1) = eps_ox; A(i,i)=-eps_ox-eps_si; A(i,i+1)=eps_si;
    elseif (i<j2)  A(i,i-1)=eps_si; A(i,i)=-2*eps_si; A(i,i+1)=1*eps_si;
    elseif (i==j2) A(i, i-1) = eps_si; A(i,i)=-eps_ox-eps_si; A(i,i+1)=eps_ox;
    elseif (i>j2) A(i,i-1)=eps_ox; A(i,i)=-2*eps_ox; A(i,i+1)=1*eps_ox;
    end
end
A(N,N)=1;    

b=zeros(N,11);
for j=1:11
    b(1,j) = 0.33374 + VG(j,1);
    for k=j1:j2
        if (k==j1+1) b(k,j) = Deltax*Deltax*q*Nacc/eps0*0.5;
    elseif (k==j2+1) b(k,j) = Deltax*Deltax*q*Nacc/eps0*0.5;
    else b(k,j)= Deltax*Deltax*q*Nacc/eps0;
        end
    end
    b(N,j) = 0.33374 + VG(j,1);
    y1(:,j)=A\b(:,j); %phi
end

x=0:Deltax:6.6e-9;

figure(1) %potential via boundary condition
plot(x, y1,'-');
xlabel('position x(m)');
ylabel('Potential (V)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');



elec= zeros(N,11);
for j=1:11
    for ii = j1:j2
    elec(ii,j) = ni*exp(q*y1(ii,j)/(k_B*T)); 
    end
end

figure(2) %electron density
plot(x, elec*1e-6,'-');

xlabel('position x(m)');
ylabel('Electron density (cm^-3)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');


c=zeros(N,11);

for j=1:11
    c(1,j) = 0.33374 + VG(j,1);
    for k=j1:j2
        if (k==j1+1) c(k,j) = Deltax*Deltax*q*(Nacc+elec(k,j))/eps0*0.5;
    elseif (k==j2+1) c(k,j) = Deltax*Deltax*q*(Nacc+elec(k,j))/eps0*0.5;
    else c(k,j)= Deltax*Deltax*q*(Nacc+elec(k,j))/eps0;
        end
    end
    c(N,j) = 0.33374 + VG(j,1);
    y3(:,j)=A\c(:,j); %phi2
end

figure(3) % updated potential phi2
plot(x, y3,'-');
xlabel('position x(m)');
ylabel('Potential (V)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');

figure(4) % difference between phi1 and phi2
plot(x, (y1-y3),'-');
xlabel('position x(m)');
ylabel('Potential (V)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');

