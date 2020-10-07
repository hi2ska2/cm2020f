clear;


nx = 9;
ny = 5;
phi = zeros(ny,nx);
N = nx*ny; 
A = zeros(N,N);
b = zeros(N,1);

for j = 1:ny
    for i = 1:nx
        n = i+nx*(j-1);
        
        if j==1     %%% Dirichlet boundary condition, bottom
            A(n,n) = 1;
            b(n,1) = 0;
            
        elseif j == ny
            if ( i <= 2)  %%% Dirichlet boundary condition, top/left
                A(n,n) = 1;
                b(n,1) = 0;
                
            elseif(i >= 8)  %%% Dirichlet boundary condition, top/right
                A(n,n) = 1;
                b(n,1) = 1;
            else  %%% Neumann boundary condition
                A(n,n-1) = 0.5;
                A(n,n+1) = 0.5;
                A(n,n) = -2;
                A(n,n-9)= 1;
                b(n,1) = 0;
            end
            
        elseif i == 1  %%% Neumann boundary condition
            A(n,n-9) = 0.5;
            A(n,n+9) = 0.5;
            A(n,n) = -2;
            A(n,n+1) = 1;
            b(n,1) = 0;
            
        elseif i == nx  %%% Neumann boundary condition
            A(n,n-9) = 0.5;
            A(n,n+9) = 0.5;
            A(n,n) = -2;
            A(n,n-1)= 1;
            b(n,1) = 0;
            
        else
            A(n,n) = -4;
            A(n,n-1) = 1;
            A(n,n+1) = 1;
            A(n,n-9) = 1;
            A(n,n+9) = 1;
            b(n,1) = 0;
        end
    end
end


x = A \ b;

for j = 1:ny
    for i = 1:nx
        phi(j,i) = x(i+nx*(j-1),1);
    end
end


surface(phi)