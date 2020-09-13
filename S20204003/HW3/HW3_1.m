L=1;   %Length of Device , m                            
N=51;
dx=L/(N-1);
A=zeros(N,N);
b=zeros(N,1);
for i=1:N
    if i==1  || i==N
        A(i,i)=1;
    else
        A(i,i-1)=1; A(i,i)=-2; A(i,i+1)=1;
    end  
end

%Boundary condition
b(1,1)=1;
b(N,1)=-1;

phi=A \ b;

x=0:dx*1e+9:L*1e+9;

%exact solution
phi_analytic=-2/L*(x/1e+9)+1;
plot(x,phi,'-o','linewidth',1.5);