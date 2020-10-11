clear;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
Deltax = 0.5e-9; % 0.5 nm spacing, dx
Deltay = 0.1e-9; % 0.1 nm spacing, dy
dydx = Deltay/Deltax;
dxdy = Deltax/Deltay;
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Nd = 1e26; % 1e20 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltax*Deltay*q/eps0;

nx = 61;  %dx = 0.5nm
ny = 71;  %dy = 0.1nm
interface1 = 11;  %y=1nm
interface2 = 61; %y=6nm

phi0 = thermal*asinh(Nd/2/ni);  %%% analytic, charge balance condition
N = nx*ny; 
X = zeros(N,1);  %tem_ElectrostaticPotential

%%% initail condition
for j = 1:ny
    for i = 1:nx
        n = i+nx*(j-1);        
        if j >= interface1 && j <= interface2
            if i < 21 || i > 41
                 X(n,1) = phi0;
            end
        end
    end
end


Vg=[0:0.02:1.1];



for v = 1:length(Vg)
    
    for newton = 1:20
        
        Jaco = sparse(N,N);
        res = zeros(N,1);
        
        for j = 1:ny
            for i = 1:nx
                n = i+nx*(j-1);
                
                if j==1
                    if i == 1
                        res(n,1) = eps_ox*(0.5*X(n+nx,1)*dxdy+0.5*X(n+1,1)*dydx-X(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx
                        res(n,1) = eps_ox*(0.5*X(n+nx,1)*dxdy+0.5*X(n-1,1)*dydx-X(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                        res(n,1) = X(n,1) - 0.33374 -Vg(v);
                        Jaco(n,n) = 1;
                        
                    else  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(X(n+nx,1)*dxdy+0.5*X(n-1,1)*dydx+0.5*X(n+1,1)*dydx-X(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j == ny
                    if i == 1
                        res(n,1) = eps_ox*(0.5*X(n-nx,1)*dxdy+0.5*X(n+1,1)*dydx-X(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx
                        res(n,1) = eps_ox*(0.5*X(n-nx,1)*dxdy+0.5*X(n-1,1)*dydx-X(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                        res(n,1) = X(n,1) - 0.33374 -Vg(v);
                        Jaco(n,n) = 1;
                        
                    else   %%% Neumann boundary condition
                        res(n,1) = eps_ox*(X(n-nx,1)*dxdy+0.5*X(n-1,1)*dydx+0.5*X(n+1,1)*dydx-X(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                    end
                    
                elseif j < interface1 || j > interface2  %%% oxide
                    if i == 1  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(0.5*X(n-nx,1)*dxdy+0.5*X(n+nx,1)*dxdy+X(n+1,1)*dydx-X(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(0.5*X(n-nx,1)*dxdy+0.5*X(n+nx,1)*dxdy+X(n-1,1)*dydx-X(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n-1) = eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    else
                        res(n,1) = eps_ox*(X(n-nx,1)*dxdy+X(n+nx,1)*dxdy+X(n-1,1)*dydx+X(n+1,1)*dydx-2*X(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -2*eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = eps_ox*dydx;
                        Jaco(n,n-1) = eps_ox*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j == interface1
                    if i == 1 || i ==nx    %%% Dirichlet boundary condition
                        res(n,1) = X(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_si*(X(n+nx,1)-X(n,1))*dxdy+eps_ox*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(-Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_si*(X(n+nx,1)-X(n,1))*dxdy+eps_ox*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(-0.5*Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    else  %%% Channel
                        res(n,1) = eps_si*(X(n+nx,1)-X(n,1))*dxdy+eps_ox*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                    end
                    
                elseif j == interface2
                    if i == 1 || i ==nx    %%% Dirichlet boundary condition
                        res(n,1) = X(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_ox*(X(n+nx,1)-X(n,1))*dxdy+eps_si*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(-Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_ox*(X(n+nx,1)-X(n,1))*dxdy+eps_si*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(-0.5*Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                        
                    else    %%% Channel
                        res(n,1) = eps_ox*(X(n+nx,1)-X(n,1))*dxdy+eps_si*(X(n-nx,1)-X(n,1))*dxdy+0.5*(eps_si+eps_ox)*(X(n-1,1)+X(n+1,1)-2*X(n,1))*dydx - 0.5*coef*(ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j > interface1 && j < interface2   %%% Silicon
                    if i == 1 || i == nx     %%% Dirichlet boundary condition, Source & Drain contact
                        res(n,1) = X(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_si*(X(n-nx,1)*dxdy+X(n+nx,1)*dxdy+X(n-1,1)*dydx+X(n+1,1)*dydx-2*X(n,1)*(dxdy+dydx)) - coef*(-Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_si*(X(n-nx,1)*dxdy+X(n+nx,1)*dxdy+X(n-1,1)*dydx+X(n+1,1)*dydx-2*X(n,1)*(dxdy+dydx)) - coef*(-0.5*Nd+ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    else   %%% Channel
                        res(n,1) = eps_si*(X(n-nx,1)*dxdy+X(n+nx,1)*dxdy+X(n-1,1)*dydx+X(n+1,1)*dydx-2*X(n,1)*(dxdy+dydx)) - coef*(ni*exp(X(n,1)/thermal)-ni*exp(-X(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(X(n,1)/thermal)+ni*exp(-X(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                    end
                    
                end
            end
        end
        
        update = Jaco \ (-res);
        X = X + update;
        
        error = max(abs(update));
        
        if error < 1e-10
            break;
        end
        
    end
end

phi = zeros(ny,nx);
for j = 1:ny
    for i = 1:nx
        phi(j,i) = X(i+nx*(j-1),1);
    end
end

x = 0:0.5:30;
y = 0:0.1:7;

surface(x,y,phi)
zlabel('Electrostatic Potential (V)')
xlabel('x (nm)')
ylabel('y (nm)')