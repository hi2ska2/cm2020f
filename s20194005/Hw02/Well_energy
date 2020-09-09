function [En] = Well_Energy(N) %input values are N=5,50,500

h_b = 6.6260 * 10^(-34)/(2*pi); % plank constant h-bar
m = 0.19 * 9.1093 * 10^(-31); % electron mass
x= 5 * 10^(-9); % length
delta_x = x/(N-1); % Length between each point
q = 1.6021* 10^(-19);

Array = []; ; % Matrix declaration

for i=2:N-3
    Array(i,i-1) = 1; 
    Array(i,i) = -2; 
    Array(i,i+1) = 1;
end

Array(1,1) = -2; 
Array(1,2) = 1;
Array(N-2,N-2) = -2;
Array(N-2,N-3) = 1; 

[Eig_V,Eig_D] = eig (Array);

d = min(Eig_D);
dm = min(abs(d));
k2 = dm/(delta_x^2); % square of k value
En = ((1/2)*(k2)*(h_b^2)/m) * 6.242*10^(18);
RealEn = ((1/2)*h_b^2)/m*(pi/(5*10^(-9)))^2/q ;
Error = (1-(En/RealEn))*100;
