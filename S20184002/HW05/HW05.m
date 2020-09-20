clear;
%*********%

% Numerical Solution
a = 0.8e-9;
b = 5e-9;
q = 1.6e-19;
Nacc = 1e25; % /m^3
N = 500;
n = round(N*(a/(2*a+b)));
delta_x = (2*a+b)/(N-1);
x = 0:delta_x:2*a+b;
xt = x';
E1 = 3.9*8.85e-12; E2 = 11.7*8.85e-12;
A = zeros(N,N);
B = zeros(N,1);
v = [1;-2;1];

for i=2:n-1
    for j=1:3
        A(i,i+j-2)=E1.*v(j,1)/(delta_x)^2;
    end
end
for i=n+1:N-n
    for j=1:3
        A(i,i+j-2)=E2.*v(j,1)/(delta_x)^2;
        B(i,1) = q*Nacc;
    end
end
for i=N-n+2:N-1
    for j=1:3
        A(i,i+j-2)=E1.*v(j,1)/(delta_x)^2;
    end
end

A(1,1) = 1; A(N,N) = 1;
A(n,n-1) = E1/(delta_x)^2; A(n,n) = -(E1+E2)/(delta_x)^2; A(n,n+1)=E2/(delta_x)^2;
A(N-n+1,N-n) = E2/(delta_x)^2; A(N-n+1,N-n+1) = -(E1+E2)/(delta_x)^2; A(N-n+1,N-n+2)=E1/(delta_x)^2;

B(n,1) = q*Nacc/2; B(N-n+1,1) = q*Nacc/2;
phi = inv(A)*B;

% A Solution
NN = 5000;
delta_xx = (2*a+b)/(NN-1);
xx = 0:delta_xx:2*a+b;
xxt = xx';
nn = round(NN*(a/(2*a+b)));
phi_N = zeros(NN,1);
for i=1:nn
    for j=1:3
        phi_N (i,1) =  -q*Nacc*b/(2*E1)*xx(1,i);
    end
end
for i=nn+1:NN-n
    for j=1:3
        phi_N (i,1) =  q*Nacc/(2*E2)*(xx(1,i)-a-b/2).^2-q*Nacc*a*b/(2*E1)-q*Nacc*b^2/(8*E2);
    end
end
for i=NN-nn+1:NN
    for j=1:3
         phi_N (i,1) =  q*Nacc*b/(2*E1)*(xx(1,i)-b-2*a);
    end
end
plot(x,phi,xx,phi_N);

