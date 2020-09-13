L=1e-3;         %Length of Device, m                      
N=501;
dx=L/(N-1);
A=zeros(N,N);
b=zeros(N,1);
middle_index=(N+1)/2;

for ii=1:N
    if ii==1  || ii==N
        A(ii,ii)=1; 
    elseif ii==middle_index
        A(ii,ii-1)=1; A(ii,ii)=-2; A(ii,ii+1)=1;
        b(ii,1)=1/(dx)*(dx)^2;
    else
        A(ii,ii-1)=1; A(ii,ii)=-2; A(ii,ii+1)=1;
    end  

end


phi=A \ b;

%exact_solution
x1=0:dx:L/2-dx;
phi_analytic1=-1/2*x1;
x2=L/2:dx:L;
phi_analytic2=1/2*x2-L/2;
phi_analytic=[transpose(phi_analytic1); transpose(phi_analytic2)];


x=0:dx*1e+9:L*1e+9;
plot(x,phi,'linewidth',1.5);