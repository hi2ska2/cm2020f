clear;
%*********%

% Numerical Solution
a = 0.8e-9;
b = 5e-9;
q = 1.602192e-19;
Nacc = 1e24; % /m^3
V_gate = 1; % gate voltage
k_B = 1.380662e-23;
T = 300;
ni = 1.075e16; % /m^3
N = 500;
n = round(N*(a/(2*a+b)));
delta_x = (2*a+b)/(N-1);
x = 0:delta_x:2*a+b;
xt = x';
E1 = 3.9*8.854187817e-12; E2 = 11.7*8.854187817e-12;
A = zeros(N,N);
B = zeros(N,1);
v = [1;-2;1];
B(N,1) = V_gate+0.33374 ;B(1,1) = V_gate+0.33374;

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
elec_con = zeros(N,1);
for i=n+1:N-n
    elec_con(i,1) = ni*exp(q*phi(i,1)/(k_B*T));
end

phi_R = zeros(N,1);
for i=2:n-1
    for j=1:3
        A(i,i+j-2)=E1.*v(j,1)/(delta_x)^2;
    end
end
for i=n+1:N-n
    for j=1:3
        A(i,i+j-2)=E2.*v(j,1)/(delta_x)^2;
        B(i,1) = q.*elec_con(i,1);
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

B(n,1) = q*elec_con(n,1)/2; B(N-n+1,1) = q*elec_con(N-n+1,1)/2;
phi_R = inv(A)*B;
diff = phi_R - phi;
plot(x,elec_con*1e-6);
figure(2)
plot(x,phi,x,phi_R);
figure(3)
plot(x,diff);


