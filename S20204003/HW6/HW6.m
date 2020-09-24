q=1.602e-19; % Elementary charge, C        
width=5e-9;         % Silicon thickness, m            
tox=0.8e-9;   % Oxide thickness, m
N=67;
dx=(width+2*tox)/(N-1);   %delta x, m
interface1=round(tox/dx)+1;  interface2=round((tox+width)/dx)+1;  
Na=1e+24;  % Doping concentration, 1/m^3
e_si=11.7;   e_ox=3.9; % Relative permittivity
e0=8.854e-12;   % Permittivity, F/m
coeff=dx*dx*q/e0;
nint=1.075e16;  % intrinsic carrier density, 1/m^3
T=300;
k_B= 1.380662e-23;  % Boltzmann constant, J/K
thermal=k_B*T/q;
eDensity=zeros(N,1);
phi=zeros(N,11);
phi_update=zeros(N,11);

for iVg=1:11
  Vg=(iVg-1)*0.1;
  A=zeros(N,N);
  b=zeros(N,1);
  for ii=1:N
   if ii==1 || ii==N
    A(ii,ii)=1;
    b(ii,1)=0.33374+Vg;
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
  %Electrostatic Potential and electron density
  phi(:,iVg)=A\b;
  eDensity(interface1:interface2,iVg)=nint*exp(phi(interface1:interface2,iVg)/thermal);


  for ii=1:N
   if ii==1 || ii==N
    A(ii,ii)=1;
    b(ii,1)=0.33374+Vg;
   elseif  interface1<ii && ii<interface2 
    A(ii,ii-1)=e_si; A(ii,ii)=-2*e_si; A(ii,ii+1)=e_si;
    b(ii,1)=coeff*(Na+eDensity(ii,iVg));
   elseif  ii<interface1 || ii>interface2
    A(ii,ii-1)=e_ox; A(ii,ii)=-2*e_ox; A(ii,ii+1)=e_ox;
    b(ii,1)=0;
   elseif ii==interface1
    A(ii,ii-1)=e_ox; A(ii,ii)=-e_ox-e_si; A(ii,ii+1)=e_si;
    b(ii,1)=coeff*(Na+eDensity(ii,iVg))/2;
   elseif ii==interface2
    A(ii,ii-1)=e_si; A(ii,ii)=-e_ox-e_si; A(ii,ii+1)=e_ox;
    b(ii,1)=coeff*(Na+eDensity(ii,iVg))/2;
   end
  end
  
   % Updated electrostatic potential
  phi_update(:,iVg)=A\b;
end

% unit change (m^-3 -> cm^-3)
eDensity=eDensity/1e6;

Difference_rate=max(abs(phi-phi_update)./phi)*100;

x=0:dx:(width+2*tox);
plot(x,phi_update,'linewidth',1.5);

