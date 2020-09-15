p1= 3.9*8.854*10^-12;
p2= 22.0*8.854*10^-12;
Thick1=0.5*10^-9;
Thick2=1.9*10^-9;

N = ((Thick1+Thick2)/10^-10)+1;
A = zeros(N,N);
for i=2:N-1;
    if i<6
        A(i,i-1)=p1; A(i,i)=-2*p1; A(i,i+1)=p1;
    end
    if i==6
         A(i,i-1)=p1; A(i,i)=-p1-p2; A(i,i+1)=p2;
    end 
    if i>6
          A(i,i-1)=p2; A(i,i)=-2*p2; A(i,i+1)=p2;  
    end
end
A(1,1)=1;
A(N,N)=1;

b=zeros(N,1);
b(N,1)=1;

x=[0:(N-1)];
y=inv(A)*b;
plot(x,y);
xlabel('x[pm]')
ylabel('V')

Cap1_analytic = p1/Thick1;
Cap2_analytic = p2/Thick2;
Cap_analytic= Cap1_analytic*Cap2_analytic/(Cap1_analytic+Cap2_analytic)*(10^-4);
Cap1_numerical = p1*(y(6,1)-y(1,1))/(0.5*10^-9)*(10^-4);
Error=(Cap1_numerical-Cap_analytic)/Cap_analytic*100;
