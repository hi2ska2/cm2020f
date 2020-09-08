q=1.602e-19;   % Elementary charge, C
h=6.626e-34; hbar=h/(2*pi); % Planck constant ,J s 
L=5e-9;       % Well length ,m
m0=9.109e-31;  m=m0*0.19;    % Electron mass (kg)  &  effective mass            
N=5;            % Discretization       
N1=N-2;          % Exclude boundary
dx=L/(N-1);        
A=zeros(N1,N1);
wavefunction=zeros(N,N1);
En_exact=zeros(N1,1);

for ii=1:N1
     if ii==1          
        A(ii,ii)=-2;  A(ii,ii+1)=1;
     elseif ii==N1             
        A(ii,ii-1)=1;  A(ii,ii)=-2;
     else                       
        A(ii,ii-1)=1; A(ii,ii)=-2; A(ii,ii+1)=1;
     end
end

[eigenvector,eigenvalue] = eig(A);
unsorted_E=-hbar*hbar*(eigenvalue)/(2*m*dx*dx);
unsorted_E=unsorted_E/q;                        % To make J to eV 
[En,ind] = sort(diag(unsorted_E));
Sorted_eigenvector=eigenvector(:,ind);


%To get normalized wavefunction
Distribution=Sorted_eigenvector.*Sorted_eigenvector;
Sum=sum(Distribution)*dx;
for ii=1:N1
wavefunction(:,ii)=[zeros(1,1) ;Sorted_eigenvector(:,ii)/sqrt(Sum(1,ii)); zeros(1,1)];  %Including boundaries
end


%Exaxt Solution of Energy and wavefunction
x=0:dx:L;
wavefunction_exact=sqrt(2/L)*sin(pi/L*x);
for n=1:N1
    En_exact(n,1)=hbar*hbar*(n*pi/L)^2/(2*m*q);
end
wavefunction_exact=transpose(wavefunction_exact);

x=linspace(0,L*1e+9,N);
plot(wavefunction(:,3),'linewidth',1.25);
xlabel('Position(nm)');
ylabel('Wave function');

