%Assignment 5, MOS
clear all;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/cm
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, cm
tsi = 5e-9; %silicon layer thickness, cm


Nacc = 1e23; %number density, m^-3
prompt = 'What is the number N?';
C1 = q*Nacc*tsi/(2*eps0*eps_ox); %analytic constant A
C2 = q*Nacc/(2*eps0*eps_si); %analytic constant B
C3 = tsi/2+tox; %analytic constant b
C4 = q*Nacc*tox*tsi/(2*eps0*eps_ox)+q*Nacc*tsi*tsi*0.5*0.5/(2*eps0*eps_si); %analytic constant C


N= input(prompt);

Deltax = (6.6e-9)/(N-1); %spacing

interface1 = round((0.8/6.6)*(N-1));
interface2 = round((5.8/6.6)*(N-1));
j1=interface1-1;
j2=interface2-1;
A=[];
A(1,1)=1;
A(N,N)=1;
y2=zeros(N,1);


for i=2:(j1);
    A(i,i-1)=eps_ox;
    A(i,i)=-2*eps_ox;
    A(i,i+1)=1*eps_ox;
    
    y2(i,1)=-Deltax*C1*(i-1); %Analytic solution
end

A(j1+1,j1)=eps_ox;
A(j1+1,j1+1)=-eps_ox-eps_si;
A(j1+1,j1+2)=eps_si;
y2(j1+1,1)=-Deltax*C1*(j1);
for i=(j1+2):(j2);
    A(i,i-1)=eps_si;
    A(i,i)=-2*eps_si;
    A(i,i+1)=1*eps_si;

     y2(i,1)=C2*(Deltax*(i-1)-C3)^2-C4; %Analytic solution
end

    A(j2+1,j2)=eps_si;
    A(j2+1,j2+1)=-eps_ox-eps_si;
    A(j2+1,j2+2)=eps_ox;
    y2(j2+1,1)=C2*(Deltax*(j2)-C3)^2-C4;
    
for i=(j2+2):(N-1);
    A(i,i-1)=eps_ox;
    A(i,i)=-2*eps_ox;
    A(i,i+1)=1*eps_ox;

    y2(i,1)=C1*(Deltax*(i-1)-2*tox-tsi);
end
    
b=zeros(N,1);
for k=j1:j2
    if (k==j1+1) b(k,1) = Deltax*Deltax*q*Nacc/eps0*0.5;
    elseif (k==j2+1) b(k,1) = Deltax*Deltax*q*Nacc/eps0*0.5;
    else b(k,1)= Deltax*Deltax*q*Nacc/eps0;
    end
end

y1=A\b;
x=0:Deltax:6.6e-9;


plot(x,y1(:,1),':r',x, y2(:,1));
xlabel('x position(m)'); 
ylabel('Potential(V)');

int1= (y2(j1,1)-y1(j1-1,1))*100/y2(j1-1,1);
midp= (y2(N/2,1)-y1(N/2,1))*100/y2(N/2,1)
int2= (y2(j2,1)-y1(j2,1))*100/y2(j2,1);
disp(sprintf ('The eerror for intersection 1 is %g %%', int1));
disp(sprintf ('The eerror for midpoint is %g %%', midp));
disp(sprintf ('The eerror for intersection 2 is %g %%', int2));