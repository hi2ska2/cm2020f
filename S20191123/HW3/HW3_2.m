clear;

N = 1001;
A = zeros(N,N);
b = zeros(N,1);

length = 1;  %length, a

dx = length/(N-1);

l = [0:1:N-1]*dx; 

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
    
    if i == (N+1)/2
        b(i,1)=dx;
    end
end



x = A \ b;


plot(l,x)

figure(2)
plot(l,b)