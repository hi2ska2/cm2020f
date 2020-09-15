%Assignment 4, capacitance
clear all;

e0 = 8.85418; % Vacuum permittivity
e1 = 3.9; % relative permittivity for material 1
e2 = 22;  % relative permittivity for material 2

prompt = 'What is the number N?';
N= input(prompt);

A=[];
A(1,1)=1;
A(N,N)=1;
j=round(5/24*N); %정수로 변환하는 함수

for i=2:(j-1);
    A(i,i-1)=e1/(1/(N-1))^2;
    A(i,i)=-2*e1/(1/(N-1))^2;
    A(i,i+1)=1*e1/(1/(N-1))^2;
    
end
 
    A(j,j-1)=e1/(1/(N-1))^2;
    A(j,j)=-e1/(1/(N-1))^2-e2/(1/(N-1))^2;
    A(j,j+1)=1*e2/(1/(N-1))^2;
for i=(j+1):(N-1);
    A(i,i-1)=e2/(1/(N-1))^2;
    A(i,i)=-2*e2/(1/(N-1))^2;
    A(i,i+1)=1*e2/(1/(N-1))^2;
    
end

b=zeros(N,1);
b(N,1)=1;


x=A\b;
y=0:24/(N-1):24;
C=e1*e0*x(j+2,1)/5;
plot(y,x(:,1));
xlabel('x position(0.1nm)'); 
ylabel('Voltage(V)');
