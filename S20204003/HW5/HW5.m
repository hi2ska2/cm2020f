q=1.602e-19; % Elementary charge, C        
width=5e-9;         % Silicon thickness, m            
tox=0.8e-9;   % Oxide thickness, m
N=67;
dx=(width+2*tox)/(N-1);   %delta x, m
interface1=round(tox/dx)+1;  interface2=round((tox+width)/dx)+1;  
Na=1e+25;  % Doping concentration, 1/m^3
e_si=11.7;   e_ox=3.9; % Relative permittivity
e0=8.854e-12;   % Permittivity, F/m
A=zeros(N,N);
b=zeros(N,1);
coeff=dx*dx*q/e0;

for ii=1:N
  if ii==1 || ii==N
    A(ii,ii)=1;
    b(ii,1)=0;
  elseif  interface1<ii && ii<interface2 
    A(ii,ii-1)=e_si; A(ii,ii)=-2*e_si; A(ii,ii+1)=e_si;
    b(ii,1)=coeff*Na;
  elseif  ii<interface1 || ii>interface2
    A(ii,ii-1)=e_ox; A(ii,ii)=-2*e_ox; A(ii,ii+1)=e_ox;
    b(ii,1)=0;
  elseif ii==interface1
    A(ii,ii-1)=e_ox; A(ii,ii)=-e_ox-e_si; A(ii,ii+1)=e_si;
    b(ii,1)=coeff*Na/2;
  elseif ii==interface2
    A(ii,ii-1)=e_si; A(ii,ii)=-e_ox-e_si; A(ii,ii+1)=e_ox;
    b(ii,1)=coeff*Na/2;
  end
    
end

phi=A\b;
x=0:dx:(width+2*tox);
plot(x,phi,'linewidth',1.5);
hold on;


%Analytic solution
phi_interface=-3*tox*q*Na*width/(2*e_si*e0);
x0=0:dx:tox;
x1=tox+dx:dx:tox+width;
x2=tox+width+dx:dx:2*tox+width;
phi0=phi_interface*x0/tox;
phi1=q*Na*(x1-tox).*(x1-tox-width)/(2*e0*e_si)+phi_interface;
phi2=-phi_interface*(x2-2*tox-width)/tox;
phi_analytic=[transpose(phi0);  transpose(phi1); transpose(phi2)];
