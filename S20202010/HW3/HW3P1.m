clear all;
prompt = 'What is the number N?';
N= input(prompt); 

A=[];
A(1,1)=1;
A(N,N)=1;

for i=2:(N-1);
    A(i,i-1)=1/(1/(N-1))^2;
    A(i,i)=-2/(1/(N-1))^2;
    A(i,i+1)=1/(1/(N-1))^2;
    
end

b=zeros(N,1);
b(N,1)=-1;
b(1,1)=1;

x=A\b;
y=0:1/(N-1):1;

plot(y,x(:,1))
xlabel('x'); 
ylabel('\Phi');


