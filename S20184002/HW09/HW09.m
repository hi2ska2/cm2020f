clear;
x = zeros(9,5);
y = zeros(9,5);
z = zeros(9,5);
phi = zeros(45,1);
A = zeros(45,45);
N = 10;

for j=1:5
    for i=1:9
        x(i,j) = i;
        y(i,j) = j;
    end
end

% Bottom position
for i=1:9
    phi(i,1) = 0;
    A(i,i) = 1;
end

% Top/Left
for i=37:38
    phi(i,1) = 0;
    A(i,i) = 1;
end

% Top/Right
for i=44:45
    phi(i,1) = 1;
    A(i,i) = 1;
end

for i=39:43
    A(i,i) = -2; A(i,i-1) = 0.5; A(i,i+1) = 0.5; A(i,i-9) = 1;
end

for i=11:17
    A(i,i) = -4; A(i,i-1) = 1; A(i,i+1) = 1; A(i,i+9) = 1; A(i,i-9) = 1;
    A(i+9,i+9) = -4; A(i+9,i+9-1) = 1; A(i+9,i+9+1) = 1; A(i+9,i+9+9) = 1; A(i+9,i+9-9) = 1;
    A(i+18,i+18) = -4; A(i+18,i+18-1) = 1; A(i+18,i+18+1) = 1; A(i+18,i+18+9) = 1; A(i+18,i+18-9) = 1;
end

for i=1:3
    A(i*9+1,i*9+1) = -2; A(i*9+1,i*9+1+1) = 1; A(i*9+1,i*9+1+9) = 0.5; A(i*9+1,i*9+1-9) = 0.5;
    A((i+1)*9,(i+1)*9) = -2; A((i+1)*9,(i+1)*9-1) = 1; A((i+1)*9,(i+1)*9+9) = 0.5; A((i+1)*9,(i+1)*9-9) = 0.5;
end

phi = A\phi;

for j=1:5
    for i=1:9
        z(i,j) = phi((j-1)*9+i,1);
    end
end

surf(x,y,z);
xlabel('X') 
ylabel('Y')
 