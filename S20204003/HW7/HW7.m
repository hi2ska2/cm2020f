q=1.602e-19; % Elementary charge, C
nint=1e16;    % Intrinsic carrier density, m-3
k_B=1.38065e-23;    % Boltzmann constant, J/K
T=300;   % Temperature, K
thermal=k_B*T/q;
phi_numerical=zeros(9,2);
phi_analytic=zeros(9,2);
update_phi=zeros(9,40);

for Case=1:2  %Case=1 :Donor Case=2 :Acceptor
   if Case==1
      polar=1;
      phi=0.3;    %initial solution
   elseif Case==2
      polar=-1;
      phi=-0.3;   %initial solution
   end
    
   for iDop=1:9
     Ndop=polar*10^(15+iDop);                % Doping concentration, m-3

     for Newton=1:40   
         res=q*(Ndop+nint*exp(-phi/thermal)-nint*exp(phi/thermal));
         Jaco=q*(-nint*exp(-phi/thermal)-nint*exp(phi/thermal))/thermal;
         update_phi =  (-res) / Jaco ;
         phi=phi+update_phi;
         
         if abs(update_phi)<1e-15
             break;
         end
     end
     phi_numerical(iDop,Case)=phi;
     %analytic solution
     phi_analytic(iDop,Case)=(thermal*asinh(Ndop/(2*nint)));
   end
end

error=abs(phi_analytic-phi_numerical)./phi_numerical*100;
