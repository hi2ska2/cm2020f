prompt = 'Choose the number of N ';
N= input(prompt);

Array = zeros(N,N); % Matrix declaration

for i=2:N-1;
    Array(i,i-1) = 1; 
    Array(i,i) = -2; 
    Array(i,i+1) = 1;
end

Array(1,1)=1;
Array(N,N)=1;

b=zeros(N,1);

if mod(N,2)==1
b((N+1)/2,1)=1;
end

x=0:1/(N-1):1;
y=inv(Array)*b;

plot(x,y);
