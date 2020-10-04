q=1.602e-19; % Elementary charge, C    
width=5e-9;    % Silicon thickness, m                   
tox=0.8e-9;   % Oxide thickness, m
k_B=1.380662e-23;  % Boltzmann constant, J/K
T=300;  % Temperature, K
N=67;    
dx=(width+2*tox)/(N-1);  
interface1=round(tox/dx)+1;  interface2=round((tox+width)/dx)+1;  
Na=1e24;   % Doping concentration, 1/m^3
nint=1.075e16; % Intrinsic carrier density, 1/m^3
e_si=11.7;  e_ox=3.9;  e0=8.854e-12;  % Permittivity, F/m
thermal=k_B*T/q; 
Elec_save=zeros(N,51);
coeff=dx*dx*q/e0;

phi=zeros(N,1);
phi(:,1)=0.33374;


for index_Vg=1:101
  Vg=(index_Vg-1)*0.01;
  for Newton=1:20
   Jaco=sparse(N,N);
   res=zeros(N,1); 
   for ii=1:N
    if ii==1 || ii==N
      res(ii,1)=phi(ii,1)-0.33374-Vg;
      Jaco(ii,ii)=1;
    elseif  interface1<ii && ii<interface2
      res(ii,1)=e_si*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal));
      Jaco(ii,ii+1)=e_si; Jaco(ii,ii)=-2*e_si-coeff*nint*exp(phi(ii,1)/thermal)/thermal;
      Jaco(ii,ii-1)=e_si;
    elseif  ii<interface1 || ii>interface2
      res(ii,1)=e_ox*(phi(ii-1,1)-2*phi(ii,1)+phi(ii+1,1));
      Jaco(ii,ii+1)=e_ox; Jaco(ii,ii)=-2*e_ox; Jaco(ii,ii-1)=e_ox;
    elseif ii==interface1
     res(ii,1)=e_si*(-phi(ii,1)+phi(ii+1,1))-e_ox*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal))/2;
     Jaco(ii,ii+1)=e_si; 
     Jaco(ii,ii)=-e_si-e_ox-coeff*nint*exp(phi(ii,1)/thermal)/(2*thermal);
     Jaco(ii,ii-1)=e_ox;
     elseif ii==interface2
     res(ii,1)=e_ox*(-phi(ii,1)+phi(ii+1,1))-e_si*(phi(ii,1)-phi(ii-1,1))+coeff*(-Na-nint*exp(phi(ii,1)/thermal))/2;
     Jaco(ii,ii+1)=e_ox; 
     Jaco(ii,ii)=-e_si-e_ox-coeff*nint*exp(phi(ii,1)/thermal)/(2*thermal);
     Jaco(ii,ii-1)=e_si;
    end
   end
   delphi=Jaco\(-res);
   phi=phi+delphi;
   
   error(Newton,index_Vg)=max(abs(delphi));
   
   if max(abs(delphi))<1e-15
     break;
   end
  end
  Elec_save(:,index_Vg)=nint*exp(phi/thermal)/1e+6;
  save_phi(:,index_Vg)=phi;   
end
Electron=[1/2*Elec_save(interface1,:);Elec_save(interface1+1:interface2-1,:);1/2*Elec_save(interface2,:)]; % Electron density, 1/cm^3

integrated_eDensity=sum(Electron)*dx*100;   % Electron density, 1/cm^2
x=0:dx:width+2*tox;
axis=0:0.01:Vg;
plot(axis,integrated_eDensity,'Linewidth',1.5);
hold on;
