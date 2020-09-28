prompt = 'What is the number N?';
N= input(prompt);
i=2;
me=0.511*10^(6)/(3.00*10^8)^2; %Electron rest mass in eV
h=6.582*10^(-16); %Planck constant (or Dirac constant), h-bar in eV * s
a=5*10^(-9); %Potential well width, (=5nm)
dx = a/(N-1); %delta x

A=[];
A(1,1)=-2;
A(1,2)=1;

if i<(N-2);
  
    A(i,i-1)=1;
    A(i,i)=-2;
    A(i,i+1)=1;
end
A(N-2,N-3) = 1;
A(N-2,N-2)= -2;

[D,e] = eig(A);
e = diag(e);

kk=-e(N-2)/(dx)^2; %k^2, k is the wavenumber
E=h*h*kk/2/me/0.19; % Unit: eV

disp(sprintf ('The energy for given N = %g is %g eV', N, E));


