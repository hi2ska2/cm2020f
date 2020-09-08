clear;

h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
m = 0.19*m0;
N = [5,10,50, 100, 300, 500];
L = 5e-9;

analytic_E = zeros(length(N),1);
E = zeros(length(N),1);
error = zeros(length(N),1);


for j = 1:length(N)    
    analytic_E(j) = hbar^2/2/m*(pi/L)^2/q; % analytic ground Energy, eV
    dx = L/(N(j)-1);
    A = zeros(N(j)-2,N(j)-2);
for i = 1:N(j)-2;
    if(i==1)
        A(i,i) =-2;
        A(i,i+1)=1;
    
    elseif(i<N(j)-2)
        A(i,i-1)=1;
        A(i,i)=-2;
        A(i,i+1)=1;
    else
        A(i,i-1)=1;
        A(i,i)=-2;
    end
end

[eigenvectors,eigenvalues] = eig(A);
a=diag(eigenvalues);
ground=max(a);

E(j) = -ground/dx^2*hbar^2/2/m/q; % numerical ground Energy, eV

P = zeros(N(j),1);
P(2:N(j)-1) = -eigenvectors(:,N(j)-2)/sqrt(sum(eigenvectors(:,N(j)-2).^2)*dx); % numerical wavefunction, m^(-1/2)
x = dx*transpose([0:N(j)-1]);
plot(x,P); hold on;

error(j) = (E(j)/analytic_E(j)-1)*100 ; % Error (numerical/ananlytic - 1) * 100, % 



end

figure(1)
analtytic_P = sqrt(2/L)*sin(pi/L*x); % analytic wavefunction, m^(-1/2)
plot(x, analtytic_P,'r');
xlabel('Position (m)')
ylabel('Wavefunction (m^{(-1/2))}')

figure(2)
plot(N,E,'o'); hold on;
plot(N,analytic_E,'r');
xlabel('N')
ylabel('Energy (eV)')

figure(3)
plot(N, error, 'o'); hold on;
xlabel('N')
ylabel('Error (%)')