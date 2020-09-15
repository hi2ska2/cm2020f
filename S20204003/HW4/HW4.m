Length1=0.5e-9;  % length for silicon oxide ,m
Length2=1.9e-9;  % length for Hafnium oxide ,m
Length=Length1+Length2;   % Full length , m
N=241;       
dx=Length/(N-1);
e1=3.9; e2=22.0;    % Relative permittivity
e0=8.854e-12;   % Vacuum permittivity, F/m
A=zeros(N,N);
interface=round(Length1/dx)+1;  %interface
Vapplied=1;
for ii=1:N
    if ii==1 || ii==N
        A(ii,ii)=1;
    elseif  ii<interface
        A(ii,ii-1)=e1;   A(ii,ii)=-2*e1;  A(ii,ii+1)=e1;
    elseif  ii>interface
        A(ii,ii-1)=e2;   A(ii,ii)=-2*e2;  A(ii,ii+1)=e2;
    elseif ii==interface
        A(ii,ii-1)=e1;   A(ii,ii)=-e2-e1; A(ii,ii+1)=e2;
    end
end
b=[zeros(N-1,1);Vapplied*ones(1,1)];
phi=A\b;

% x=0:dx:Length;
% plot(x,phi,'linewidth',2);

DP_A=zeros(N,N);

%Matrix for displacement field

for ii=1:N
    if ii==1
        DP_A(ii,ii)=-e1*e0/dx;  DP_A(ii,ii+1)=e1*e0/dx;   
    elseif ii==N
        DP_A(ii,ii)=e2*e0/dx;  DP_A(ii,ii-1)=-e2*e0/dx;
    elseif  ii<interface
        DP_A(ii,ii)=-e1*e0/dx;  DP_A(ii,ii+1)=e1*e0/dx; 
    elseif  ii>interface
        DP_A(ii,ii)=e2*e0/dx;  DP_A(ii,ii-1)=-e2*e0/dx;
    elseif ii==interface
        DP_A(ii,ii)=e1*e0/dx;  DP_A(ii,ii-1)=-e1*e0/dx; 
    end
end


%   /1e+4 is for unit change(m^2 ->cm^2) 
Displacement=(-DP_A*phi)/Vapplied/1e+4;
C=-Displacement(1,1);


% Analytical way
% For Potential
x=0:dx:Length;
x1=0:dx:Length1;
x2=Length1+dx:dx:Length;
phi1_analytic=Vapplied/(Length1+Length2*(e1/e2))*x1;
phi2_analytic=Vapplied/((e2/e1)*Length1+Length2)*(x2-Length)+Vapplied; 
phi_analytic=[transpose(phi1_analytic); transpose(phi2_analytic)];

% For Capacitance
%   /1e+4 is for unit change(m^2 ->cm^2) 
analytical_c1=e1*e0/Length1/1e+4;
analytical_c2=e2*e0/Length2/1e+4;
analytical_C=1/(1/analytical_c1+1/analytical_c2);

Error=(analytical_C-C)/analytical_C*100;
