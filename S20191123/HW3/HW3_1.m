clear;

N = 5;
A = zeros(N,N);

for i = 1:N
    if i==1
        A(1,1) = 1;
    elseif i==N
        A(N,N)=1;
    else
        A(i,i) = -2;
        A(i,i-1) = 1;
        A(i,i+1) = 1;
    end
end

b = zeros(N,1);
b(1,1)=1;
b(N,1)=-1;

x = A \ b;

l = [0:1:N-1];
plot(l,x)