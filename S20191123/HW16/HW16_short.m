clear;


for r = 1:3
    if r == 1 % 5nm spacing
        Deltax = 5e-9;
        N = 13; % 60-nm-long structure
        interface1 = 3; % At x=10 nm
        interface2 = 11; % At x=50 nm
        
    elseif r == 2 % 1nm spacing
        Deltax = 1e-9; % 1 nm spacing
        N = 61; % 60-nm-long structure
        interface1 = 11; % At x=10 nm
        interface2 = 51; % At x=50 nm
        
    else % 0.2nm spacing
        Deltax = 0.2e-9; % 0.2 nm spacing
        N = 301; % 60-nm-long structure
        interface1 = 51; % At x=10 nm
        interface2 = 251; % At x=50 nm
    end
    
    q = 1.602192e-19; % Elementary charge, C
    eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
    k_B = 1.380662e-23; % Boltzmann constant, J/K
    T = 300.0; % Temperature, K
    thermal = k_B*T/q; % Thermal voltage, V
    x = Deltax*transpose([0:N-1]); % real space, m
    eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
    Ndon = 2e23*ones(N,1); % 2e17 /cm^3
    Ndon(1:interface1,1) = 5e25; % 5e19 /cm^3
    Ndon(interface2:N,1) = 5e25; % 5e19 /cm^3
    ni = 1.075e16; % 1.075e10 /cm^3
    coef = Deltax*Deltax*q/eps0;
    mu = 1000e-4; %mobility, 1000cm^2/vs
    A = 1e-12; %Area, 1um^2
    
    
    phi = zeros(N,1);
    phi(:,1) = thermal*log(Ndon(:,1)/ni);
    for newton=1:10
        res = zeros(N,1);
        Jaco = sparse(N,N);
        res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
        Jaco(1,1) = 1.0;
        for ii=2:N-1
            res(ii,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1));
            Jaco(ii,ii-1) = eps_si;
            Jaco(ii,ii ) = -2*eps_si;
            Jaco(ii,ii+1) = eps_si;
        end
        res(N,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
        Jaco(N,N) = 1.0;
        for ii=2:N-1
            res(ii,1) = res(ii,1) - coef*(-Ndon(ii,1)+ni*exp(phi(ii,1)/thermal));
            Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,1)/thermal)/thermal;
        end
        update = Jaco \ (-res);
        phi = phi + update;
        norm(update,inf)
    end
    
    
    
    %%% continuity equation  %%%
    elec = zeros(N,1);
    res_elec = zeros(N,1);
    Jaco_elec = sparse(N,N);
    
    res_elec(1,1) = elec(1,1) - Ndon(1,1);
    Jaco_elec(1,:) = 0.0;
    Jaco_elec(1,1) = 1.0;
    
    for ii=2:N-1
        res_elec(ii,1) = [0.5*(elec(ii+1,1)+elec(ii,1)) * (phi(ii+1,1)-phi(ii,1))/Deltax - thermal * (elec(ii+1,1)-elec(ii,1))/Deltax] - [0.5*(elec(ii,1)+elec(ii-1,1)) * (phi(ii,1)-phi(ii-1,1))/Deltax - thermal * (elec(ii,1)-elec(ii-1,1))/Deltax];
        Jaco_elec(ii,ii+1) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax - thermal / Deltax;
        Jaco_elec(ii,ii ) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax + thermal / Deltax - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax + thermal / Deltax;
        Jaco_elec(ii,ii-1) = - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax - thermal / Deltax;
    end
    
    
    res_elec(N,1) = elec(N,1) - Ndon(N,1);
    Jaco_elec(N,:) = 0.0;
    Jaco_elec(N,N) = 1.0;
    
    update_elec = Jaco_elec \ (-res_elec);
    elec = elec + update_elec;
    
    
    %%% I-V %%%
    
    for bias=1:11
        V_applied(bias,r) = 0.05 * (bias-1);
        
        for newton=1:10
            res = zeros(2*N,1);
            Jaco = sparse(2*N,2*N);
            
            %%% Poisson %%%
            res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
            Jaco(1,1) = 1.0;
            
            for ii=2:N-1
                res(2*ii-1,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1)) + coef*(Ndon(ii,1)-elec(ii,1));
                Jaco(2*ii-1,2*ii+1) = eps_si;
                Jaco(2*ii-1,2*ii-1) = -2*eps_si;
                Jaco(2*ii-1,2*ii-3) = eps_si;
                Jaco(2*ii-1,2*ii ) = -coef;
            end
            res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni) - V_applied(bias,r);
            Jaco(2*N-1,2*N-1) = 1.0;
            
            %%% Continuity %%%
            for ii=2:N-1 % edge-wise construction
                res(2*ii,1) = [0.5*(elec(ii+1,1)+elec(ii,1)) * (phi(ii+1,1)-phi(ii,1))/Deltax - thermal * (elec(ii+1,1)-elec(ii,1))/Deltax] - [0.5*(elec(ii,1)+elec(ii-1,1)) * (phi(ii,1)-phi(ii-1,1))/Deltax - thermal * (elec(ii,1)-elec(ii-1,1))/Deltax];
                Jaco(2*ii,2*ii+2) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax - thermal / Deltax;
                Jaco(2*ii,2*ii ) = 0.5* (phi(ii+1,1)-phi(ii,1))/Deltax + thermal / Deltax - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax + thermal / Deltax;
                Jaco(2*ii,2*ii-2) = - 0.5* (phi(ii,1)-phi(ii-1,1))/Deltax - thermal / Deltax;
                Jaco(2*ii,2*ii+1) = [0.5*(elec(ii+1,1)+elec(ii,1)) /Deltax];
                Jaco(2*ii,2*ii-1) = [ - 0.5*(elec(ii+1,1)+elec(ii,1))/Deltax] - [0.5*(elec(ii,1)+elec(ii-1,1)) /Deltax];
                Jaco(2*ii,2*ii-3) = [0.5*(elec(ii,1)+elec(ii-1,1)) /Deltax];
                
                
            end
            
            res(2,1) = elec(1,1) - Ndon(1,1);
            Jaco(2,:) = 0.0;
            Jaco(2,2) = 1.0;
            res(2*N,1) = elec(N,1) - Ndon(N,1);
            Jaco(2*N,:) = 0.0;
            Jaco(2*N,2*N) = 1.0;
            
            Cvector = zeros(2*N,1);
            Cvector(1:2:2*N-1,1) = thermal;
            Cvector(2:2:2*N ,1) = max(abs(Ndon));
            Cmatrix = spdiags(Cvector,0,2*N,2*N);
            Jaco_scaled = Jaco * Cmatrix;
            Rvector = 1./sum(abs(Jaco_scaled),2);
            Rmatrix = spdiags(Rvector,0,2*N,2*N);
            Jaco_scaled = Rmatrix * Jaco_scaled;
            res_scaled = Rmatrix * res;
            update_scaled = Jaco_scaled \ (-res_scaled);
            update = Cmatrix * update_scaled;
            
            phi = phi + update(1:2:2*N-1,1);
            elec = elec + update(2:2:2*N,1);
            norm(update(1:2:2*N-1,1),inf)
        end
        
        I(bias,r) = q*mu*(elec(N,1)*(phi(N,1)-phi(N-1,1))/Deltax-thermal*(elec(N,1)-elec(N-1,1))/Deltax) * A;
        
    end
    
    if r == 1
        plot(V_applied(:,1),I(:,1),'r');
        hold on;
    elseif r == 2
        plot(V_applied(:,1),I(:,2),'b');
        hold on;
    else
        plot(V_applied(:,1),I(:,3),'g');
        hold on;
    end
    
end

xlabel('Bias Voltage (V)')
ylabel('I (A)')

