clear;
%*********%
N = 10;
a = 1;
delta_x = a/(N-1);
x = 0:delta_x:1; 

A = zeros(N,N);
v = [1;-2;1];
for i=2:N-1
    for j=1:3
        A(i,i+j-2)=v(j,1)/(delta_x)^2;
    end
end
A(1,1)=1; A(N,N)=1;

B = zeros(N,1);
B(1,1) = 1; B(N,1) = -1;

phi = inv(A)*B;
plot(x,phi);
