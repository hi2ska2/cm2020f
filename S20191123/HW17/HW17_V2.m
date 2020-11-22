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
mu = 1417e-4; %mobility, 1417cm^2/vs
W = 1e-6; %Area, 1um

nx = 61;  %dx = 0.5nm
ny = 71;  %dy = 0.1nm
interface1 = 11;  %y=1nm
interface2 = 61; %y=6nm

VG = 1.1; % Gate voltage
VD = 1; % Drain voltage

phi0 = thermal*asinh(Nd/2/ni);  %%% analytic, charge balance  condition
N = nx*ny;
phi = zeros(N,1);  %ElectrostaticPotential
phi_temp = zeros(N,1);  %ElectrostaticPotential
elec = zeros(N,1);
hole = zeros(N,1);

%%% initial condition
for j = 1:ny
    for i = 1:nx
        n = i+nx*(j-1);
        if j >= interface1 && j <= interface2
            if i < 21 || i > 41
                phi(n,1) = phi0;
            end
        end
    end
end



Vg=[0:0.1:VG];
Vd=transpose([0:0.05:VD]); % 0.05V step size
%%% Poisson %%%


for v = 1:length(Vg)
    
    if v > 1
        phi = phi_temp; % 이전 step에서의 potential
    end
    
    for newton = 1:20
        
        Jaco = sparse(N,N);
        res = zeros(N,1);
        
        for j = 1:ny
            for i = 1:nx
                n = i+nx*(j-1);
                
                if j==1
                    if i == 1
                        res(n,1) = eps_ox*(0.5*phi(n+nx,1)*dxdy+0.5*phi(n+1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx
                        res(n,1) = eps_ox*(0.5*phi(n+nx,1)*dxdy+0.5*phi(n-1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                        res(n,1) = phi(n,1) - 0.33374 -Vg(v);
                        Jaco(n,n) = 1;
                        
                    else  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(phi(n+nx,1)*dxdy+0.5*phi(n-1,1)*dydx+0.5*phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j == ny
                    if i == 1
                        res(n,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx
                        res(n,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n-1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*0.5*(dxdy+dydx);
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                        res(n,1) = phi(n,1) - 0.33374 -Vg(v);
                        Jaco(n,n) = 1;
                        
                    else   %%% Neumann boundary condition
                        res(n,1) = eps_ox*(phi(n-nx,1)*dxdy+0.5*phi(n-1,1)*dydx+0.5*phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-1) = 0.5*eps_ox*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                    end
                    
                elseif j < interface1 || j > interface2  %%% oxide
                    if i == 1  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+nx,1)*dxdy+phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    elseif i == nx  %%% Neumann boundary condition
                        res(n,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+nx,1)*dxdy+phi(n-1,1)*dydx-phi(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -eps_ox*(dxdy+dydx);
                        Jaco(n,n-1) = eps_ox*dydx;
                        Jaco(n,n-nx) = 0.5*eps_ox*dxdy;
                        Jaco(n,n+nx) = 0.5*eps_ox*dxdy;
                        
                    else
                        res(n,1) = eps_ox*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx));
                        Jaco(n,n) = -2*eps_ox*(dxdy+dydx);
                        Jaco(n,n+1) = eps_ox*dydx;
                        Jaco(n,n-1) = eps_ox*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j == interface1
                    if i == 1 || i ==nx    %%% Dirichlet boundary condition
                        res(n,1) = phi(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-0.5*Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    else  %%% Channel
                        res(n,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_ox*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                    end
                    
                elseif j == interface2
                    if i == 1 || i ==nx    %%% Dirichlet boundary condition
                        res(n,1) = phi(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-0.5*Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                        
                    else    %%% Channel
                        res(n,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -(eps_si+eps_ox)*(dxdy+dydx) - 0.5*coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-1) = 0.5*(eps_si+eps_ox)*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_ox*dxdy;
                    end
                    
                elseif j > interface1 && j < interface2   %%% Silicon
                    if i == 1 || i == nx     %%% Dirichlet boundary condition, Source & Drain contact
                        res(n,1) = phi(n,1) - phi0;
                        Jaco(n,n) = 1;
                        
                    elseif i < 21 || i > 41    %%%Souce, Drain
                        res(n,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(-Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    elseif i == 21 || i == 41    %%%Souce, Drain interface
                        res(n,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(-0.5*Nd+ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                        
                    else   %%% Channel
                        res(n,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(ni*exp(phi(n,1)/thermal)-ni*exp(-phi(n,1)/thermal));
                        Jaco(n,n) = -2*eps_si*(dxdy+dydx) - coef*(ni*exp(phi(n,1)/thermal)+ni*exp(-phi(n,1)/thermal))/thermal;
                        Jaco(n,n+1) = eps_si*dydx;
                        Jaco(n,n-1) = eps_si*dydx;
                        Jaco(n,n-nx) = eps_si*dxdy;
                        Jaco(n,n+nx) = eps_si*dxdy;
                    end
                    
                end
            end
        end
        
        update = Jaco \ (-res);
        phi = phi + update;
        
        error = max(abs(update));
        
        if error < 1e-10
            break;
        end
        
    end
    
    
    phi_temp = phi;  % equilibrium에서의 potential
    
    %%% initial electron
    for j = 1:ny
        for i = 1:nx
            n = i+nx*(j-1);
            if j >= interface1 && j <= interface2
                elec(n,1) = ni*exp(phi(n,1)/thermal);
            end
        end
    end
    
    %%% initial hole
    for j = 1:ny
        for i = 1:nx
            n = i+nx*(j-1);
            if j >= interface1 && j <= interface2
                hole(n,1) = ni*exp(-phi(n,1)/thermal);
            end
        end
    end
    
    
    %%% Coupling Poisson-Scharffeter Gummel
    
    I=zeros(length(Vd),1); % Drain current
    
    for bias=1:length(Vd)
        %     V_applied(bias,1) = 0.05 * (bias-1);
        
        for newton = 1:10
            
            Jaco = sparse(3*N,3*N);
            res = zeros(3*N,1);
            
            %%% Poisson %%%
            
            for j = 1:ny
                for i = 1:nx
                    n = i+nx*(j-1);
                    
                    if j==1
                        if i == 1
                            res(3*n-2,1) = eps_ox*(0.5*phi(n+nx,1)*dxdy+0.5*phi(n+1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*0.5*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2+nx*3) = 0.5*eps_ox*dxdy;
                            
                        elseif i == nx
                            res(3*n-2,1) = eps_ox*(0.5*phi(n+nx,1)*dxdy+0.5*phi(n-1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*0.5*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2-3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2+nx*3) = 0.5*eps_ox*dxdy;
                            
                        elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                            res(3*n-2,1) = phi(n,1) - 0.33374 -Vg(v);
                            Jaco(3*n-2,3*n-2) = 1;
                            
                        else  %%% Neumann boundary condition
                            res(3*n-2,1) = eps_ox*(phi(n+nx,1)*dxdy+0.5*phi(n-1,1)*dydx+0.5*phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_ox*dxdy;
                        end
                        
                    elseif j == ny
                        if i == 1
                            res(3*n-2,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*0.5*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = 0.5*eps_ox*dxdy;
                            
                        elseif i == nx
                            res(3*n-2,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n-1,1)*dydx-phi(n,1)*0.5*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*0.5*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2-3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = 0.5*eps_ox*dxdy;
                            
                        elseif i>=21 && i<=41  %%% gate, Dirichlet boundary condition
                            res(3*n-2,1) = phi(n,1) - 0.33374 -Vg(v);
                            Jaco(3*n-2,3*n-2) = 1;
                            
                        else   %%% Neumann boundary condition
                            res(3*n-2,1) = eps_ox*(phi(n-nx,1)*dxdy+0.5*phi(n-1,1)*dydx+0.5*phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_ox*dxdy;
                        end
                        
                    elseif j < interface1 || j > interface2  %%% oxide
                        if i == 1  %%% Neumann boundary condition
                            res(3*n-2,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+nx,1)*dxdy+phi(n+1,1)*dydx-phi(n,1)*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = 0.5*eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = 0.5*eps_ox*dxdy;
                            
                        elseif i == nx  %%% Neumann boundary condition
                            res(3*n-2,1) = eps_ox*(0.5*phi(n-nx,1)*dxdy+0.5*phi(n+nx,1)*dxdy+phi(n-1,1)*dydx-phi(n,1)*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -eps_ox*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2-3) = eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = 0.5*eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = 0.5*eps_ox*dxdy;
                            
                        else
                            res(3*n-2,1) = eps_ox*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx));
                            Jaco(3*n-2,3*n-2) = -2*eps_ox*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-3) = eps_ox*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_ox*dxdy;
                        end
                        
                    elseif j == interface1
                        if i == 1     %%% Source contact
                            res(3*n-2,1) = phi(n,1) - phi0;
                            Jaco(3*n-2,3*n-2) = 1;
                            
                            
                        elseif i ==nx    %%% Drain contact
                            res(3*n-2,1) = phi(n,1) - phi0 -Vd(bias,1);
                            Jaco(3*n-2,3*n-2) = 1;
                            
                        elseif i < 21 || i > 41    %%%Souce, Drain
                            res(3*n-2,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                            
                        elseif i == 21 || i == 41    %%%Souce, Drain interface
                            res(3*n-2,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-0.5*Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                            
                        else  %%% Channel
                            res(3*n-2,1) = eps_si*(phi(n+nx,1)-phi(n,1))*dxdy+eps_ox*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                        end
                        
                    elseif j == interface2
                        if i == 1     %%% Source contact
                            res(3*n-2,1) = phi(n,1) - phi0;
                            Jaco(3*n-2,3*n-2) = 1;
                            
                            
                        elseif i ==nx    %%% Drain contact
                            res(3*n-2,1) = phi(n,1) - phi0 -Vd(bias,1);
                            Jaco(3*n-2,3*n-2) = 1;
                            
                        elseif i < 21 || i > 41    %%%Souce, Drain
                            res(3*n-2,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                            
                        elseif i == 21 || i == 41    %%%Souce, Drain interface
                            res(3*n-2,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(-0.5*Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                            
                        else    %%% Channel
                            res(3*n-2,1) = eps_ox*(phi(n+nx,1)-phi(n,1))*dxdy+eps_si*(phi(n-nx,1)-phi(n,1))*dxdy+0.5*(eps_si+eps_ox)*(phi(n-1,1)+phi(n+1,1)-2*phi(n,1))*dydx - 0.5*coef*(elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -(eps_si+eps_ox)*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-3) = 0.5*(eps_si+eps_ox)*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_ox*dxdy;
                            Jaco(3*n-2,3*n-1) = - 0.5*coef;
                            Jaco(3*n-2,3*n) =  0.5*coef;
                        end
                        
                    elseif j > interface1 && j < interface2   %%% Silicon
                        if i == 1     %%% Source contact
                            res(3*n-2,1) = phi(n,1) - phi0;
                            Jaco(3*n-2,3*n-2) = 1;
                            
                            
                        elseif i ==nx    %%% Drain contact
                            res(3*n-2,1) = phi(n,1) - phi0 -Vd(bias,1);
                            Jaco(3*n-2,3*n-2) = 1;
                            
                        elseif i < 21 || i > 41    %%%Souce, Drain
                            res(3*n-2,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(-Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -2*eps_si*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - coef;
                            Jaco(3*n-2,3*n) =  coef;
                            
                        elseif i == 21 || i == 41    %%%Souce, Drain interface
                            res(3*n-2,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(-0.5*Nd+elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -2*eps_si*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - coef;
                            Jaco(3*n-2,3*n) =  coef;
                            
                        else   %%% Channel
                            res(3*n-2,1) = eps_si*(phi(n-nx,1)*dxdy+phi(n+nx,1)*dxdy+phi(n-1,1)*dydx+phi(n+1,1)*dydx-2*phi(n,1)*(dxdy+dydx)) - coef*(elec(n,1)-hole(n,1));
                            Jaco(3*n-2,3*n-2) = -2*eps_si*(dxdy+dydx);
                            Jaco(3*n-2,3*n-2+3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-3) = eps_si*dydx;
                            Jaco(3*n-2,3*n-2-nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-2+nx*3) = eps_si*dxdy;
                            Jaco(3*n-2,3*n-1) = - coef;
                            Jaco(3*n-2,3*n) =  coef;
                        end
                        
                    end
                end
            end
            
            
            %%% Scharffeter Gummel %%%
            for j = interface1:interface2
                for i = 1:nx
                    n = i+nx*(j-1);
                    
                    if i >=2 && i<=nx-1
                        x1=(phi(n+1,1)-phi(n,1))/thermal;
                        x2=(phi(n,1)-phi(n+1,1))/thermal;
                        x3=(phi(n,1)-phi(n-1,1))/thermal;
                        x4=(phi(n-1,1)-phi(n,1))/thermal;
                        
                        if abs((phi(n+1,1)-phi(n,1))/thermal)< 2.502000000000000e-2 || abs((phi(n,1)-phi(n-1,1))/thermal)< 2.502000000000000e-2
                            BP1 = 1 - x1/2 + x1^2/12*(1-x1^2/60*(1-x1^2/42));
                            BN1 = 1 - x2/2 + x2^2/12*(1-x2^2/60*(1-x2^2/42));
                            BP2 = 1 - x3/2 + x3^2/12*(1-x3^2/60*(1-x3^2/42));
                            BN2 = 1 - x4/2 + x4^2/12*(1-x4^2/60*(1-x4^2/42));
                            
                            dBP1 = -0.5 + x1/6*(1-x1^2/30*(1-x1^2/28));
                            dBN1 = -0.5 + x2/6*(1-x2^2/30*(1-x2^2/28));
                            dBP2 = -0.5 + x3/6*(1-x3^2/30*(1-x3^2/28));
                            dBN2 = -0.5 + x4/6*(1-x4^2/30*(1-x4^2/28));
                            
                        elseif abs((phi(n+1,1)-phi(n,1))/thermal)< 1.500000000000000e-1 || abs((phi(n,1)-phi(n-1,1))/thermal)<  1.500000000000000e-1
                            BP1 = 1 - x1/2 + x1^2/12*(1-x1^2/60*(1-x1^2/42*(1-x1^2/40*(1-0.02525252525252525252525*x1^2))));
                            BN1 = 1 - x2/2 + x2^2/12*(1-x2^2/60*(1-x2^2/42*(1-x2^2/40*(1-0.02525252525252525252525*x2^2))));
                            BP2 = 1 - x3/2 + x3^2/12*(1-x3^2/60*(1-x3^2/42*(1-x3^2/40*(1-0.02525252525252525252525*x3^2))));
                            BN2 = 1 - x4/2 + x4^2/12*(1-x4^2/60*(1-x4^2/42*(1-x4^2/40*(1-0.02525252525252525252525*x4^2))));
                            
                            dBP1 = -0.5 + x1/6*(1-x1^2/30*(1-x1^2/28*(1-x1^2/30*(1-0.03156565656565656565657*x1^2))));
                            dBN1 = -0.5 + x2/6*(1-x2^2/30*(1-x2^2/28*(1-x2^2/30*(1-0.03156565656565656565657*x2^2))));
                            dBP2 = -0.5 + x3/6*(1-x3^2/30*(1-x3^2/28*(1-x3^2/30*(1-0.03156565656565656565657*x3^2))));
                            dBN2 = -0.5 + x4/6*(1-x4^2/30*(1-x4^2/28*(1-x4^2/30*(1-0.03156565656565656565657*x4^2))));
                            
                        else
                            BP1 = x1/(exp(x1)-1);
                            BN1 = x2/(exp(x2)-1);
                            BP2 = x3/(exp(x3)-1);
                            BN2 = x4/(exp(x4)-1);
                            
                            dBP1 = 1/(exp(x1)-1)-x1*exp(x1)/(exp(x1)-1)^2;
                            dBN1 = 1/(exp(x2)-1)-x2*exp(x2)/(exp(x2)-1)^2;
                            dBP2 = 1/(exp(x3)-1)-x3*exp(x3)/(exp(x3)-1)^2;
                            dBN2 = 1/(exp(x4)-1)-x4*exp(x4)/(exp(x4)-1)^2;
                            
                        end
                        
                        %%% electron %%%
                        
                        res(3*n-1,1) = elec(n+1,1)*BP1 - elec(n,1)*BN1 - elec(n,1)*BP2 + elec(n-1,1)*BN2;
                        Jaco(3*n-1,3*n-1+3) = BP1;
                        Jaco(3*n-1,3*n-1 ) = -BN1 -BP2;
                        Jaco(3*n-1,3*n-1-3) = BN2;
                        Jaco(3*n-1,3*n-1+2) = elec(n+1,1)*dBP1/thermal + elec(n,1)*dBN1/thermal;
                        Jaco(3*n-1,3*n-1-1) = -elec(n+1,1)*dBP1/thermal -elec(n,1)*dBN1/thermal -elec(n,1)*dBP2/thermal -elec(n-1,1)*dBN2/thermal;
                        Jaco(3*n-1,3*n-1-4) = elec(n,1)*dBP2/thermal + elec(n-1,1)*dBN2/thermal;
                        
                        
                        %%% hole %%%
                        
                        res(3*n,1) = hole(n+1,1)*BN1 - hole(n,1)*BP1 - hole(n,1)*BN2 + hole(n-1,1)*BP2;
                        Jaco(3*n,3*n+3) = BN1;
                        Jaco(3*n,3*n ) = -BP1 -BN2;
                        Jaco(3*n,3*n-3) = BP2;
                        Jaco(3*n,3*n+1) = -hole(n+1,1)*dBN1/thermal - hole(n,1)*dBP1/thermal;
                        Jaco(3*n,3*n-2) = hole(n+1,1)*dBN1/thermal +hole(n,1)*dBP1/thermal +hole(n,1)*dBN2/thermal +hole(n-1,1)*dBP2/thermal;
                        Jaco(3*n,3*n-5) = -hole(n,1)*dBN2/thermal - hole(n-1,1)*dBP2/thermal;
                        
                        %             elseif i == 1
                        %                 res(3*n-1,1) = elec(n,1) - Nd;
                        %                 Jaco(3*n-1,:) = 0.0;
                        %                 Jaco(3*n-1,3*n-1) = 1.0;
                        
                    else
                        
                        %%% electron %%%
                        
                        res(3*n-1,1) = elec(n,1) - Nd;
                        Jaco(3*n-1,:) = 0.0;
                        Jaco(3*n-1,3*n-1) = 1.0;
                        
                        %%% hole %%%
                        
                        res(3*n,1) = hole(n,1) - ni^2/Nd;
                        Jaco(3*n,:) = 0.0;
                        Jaco(3*n,3*n) = 1.0;
                        
                    end
                    
                end
                
            end
            
            
            iii = find(sum(abs(Jaco),2)==0);  % electron과 hole Jacobian에서 Oxide부분의 0이 되는 부분들을 보정.(NaN이 안 나오게)
            Jaco_helper = sparse(iii,iii,1,3*N,3*N);
            Jaco = Jaco + Jaco_helper;
            
            Cvector = zeros(3*N,1);
            Cvector(1:3:3*N-2,1) = thermal;
            Cvector(2:3:3*N-1 ,1) = max(abs(Nd));
            Cvector(3:3:3*N ,1) = max(abs(ni^2/Nd));
            Cmatrix = spdiags(Cvector,0,3*N,3*N);
            Jaco_scaled = Jaco * Cmatrix;
            Rvector = 1./sum(abs(Jaco_scaled),2);
            Rmatrix = spdiags(Rvector,0,3*N,3*N);
            Jaco_scaled = Rmatrix * Jaco_scaled;
            res_scaled = Rmatrix * res;
            update_scaled = Jaco_scaled \ (-res_scaled);
            update = Cmatrix * update_scaled;
            
            
            phi = phi + update(1:3:3*N-2,1);
            elec = elec + update(2:3:3*N-1,1);
            hole = hole + update(3:3:3*N,1);
            %         norm(update(1:3:3*N-2,1),inf)
            
            if norm(update(1:3:3*N-2,1),inf) < 1e-10
                break;
            end
            
        end
        
        for j = interface1:interface2
            i = nx;
            n = i+nx*(j-1);
            if j == interface1 || j == interface2
                I(bias,1) = I(bias,1) + q*mu*(elec(n,1)*(phi(n,1)-phi(n-1,1))/Deltax-thermal*(elec(n,1)-elec(n-1,1))/Deltax) * W*Deltay/2;
                
            else
                I(bias,1) = I(bias,1) + q*mu*(elec(n,1)*(phi(n,1)-phi(n-1,1))/Deltax-thermal*(elec(n,1)-elec(n-1,1))/Deltax) * W*Deltay;
                
            end
            
        end
    end
    
    
    plot(Vd,I); hold on;
end



xlabel('Drain Voltage (V)')
ylabel('I (A)')

Phi = zeros(ny,nx); % electrostatic potential 2D
for j = 1:ny
    for i = 1:nx
        Phi(j,i) = phi(i+nx*(j-1),1);
    end
end


Elec = zeros(ny,nx); % electron 2D
for j = 1:ny
    for i = 1:nx
        Elec(j,i) = elec(i+nx*(j-1),1);
    end
end

Hole = zeros(ny,nx); % hole 2D
for j = 1:ny
    for i = 1:nx
        Hole(j,i) = hole(i+nx*(j-1),1);
    end
end

% x = 0:0.5:30;
% y = 0:0.1:7;
%
% surface(x,y,Phi)
% zlabel('Electrostatic Potential (V)')
% xlabel('x (nm)')
% ylabel('y (nm)')