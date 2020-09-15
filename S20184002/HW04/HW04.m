clear;
%*********%

% Numerical Solution
a = 0.5e-9;
b = 1.9e-9;
N = 500;
n = round(N*(a/(a+b)));
delta_x = (a+b)/(N-1);
x = 0:delta_x:a+b;
E1 = 3.9*8.85e-12; E2 = 22*8.85e-12;
A = zeros(N,N);
v = [1;-2;1];

for i=2:n-1
    for j=1:3
        A(i,i+j-2)=E1.*v(j,1)/(delta_x)^2;
    end
end
for i=n+1:N-1
    for j=1:3
        A(i,i+j-2)=E2.*v(j,1)/(delta_x)^2;
    end
end

A(1,1) = 1; A(N,N) = 1;
A(n,n-1) = E1/(delta_x)^2; A(n,n) = -(E1+E2)/(delta_x)^2; A(n,n+1)=E2/(delta_x)^2;
B = zeros(N,1); B(N,1) = 1;
phi = inv(A)*B;

% Analytic solution 
NN = 5000;
nn = round(NN*(a/(a+b)));
xx = 0:(a+b)/(NN-1):a+b;
An = zeros(NN,1); 
c = E2/(E1*b+E2*a);
d = E1/(E1*b+E2*a);
e = (E2*a-E1*a)/(E1*b+E2*a);

for i=1:nn
    An(i,1) = c*xx(1,i);
end
for i=nn:NN
    An(i,1) = d*xx(1,i)+e;
end

% Analytic capacitance
C1 = E1/a; C2 = E2/b; C = 1/((1/C1)+(1/C2)); 

plot(x,phi,xx,An);

