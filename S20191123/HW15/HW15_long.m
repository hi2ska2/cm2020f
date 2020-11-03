clear;


for r = 1:3
    if r == 1 % 10nm spacing
        Deltax = 10e-9; % 1 nm spacing // 0.5nm spacing  // 10nm spacing
        N = 61; % 600-nm-long structure //601  // 1201  //  61
        interface1 = 11; % At x=100 nm //101 // 201  // 11
        interface2 = 51; % At x=500 nm //501 // 1001 // 51
        
    elseif r == 2 % 1nm spacing
        Deltax = 1e-9; % 1 nm spacing // 0.5nm spacing  // 10nm spacing
        N = 601; % 600-nm-long structure //601  // 1201  //  61
        interface1 = 101; % At x=100 nm //101 // 201  // 11
        interface2 = 501; % At x=500 nm //501 // 1001 // 51
        
    else % 0.5nm spacing
        Deltax = 0.5e-9; % 1 nm spacing // 0.5nm spacing  // 10nm spacing
        N = 1201; % 600-nm-long structure //601  // 1201  //  61
        interface1 = 201; % At x=100 nm //101 // 201  // 11
        interface2 = 1001; % At x=500 nm //501 // 1001 // 51
    end
    
    q = 1.602192e-19; % Elementary charge, C
    eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
    k_B = 1.380662e-23; % Boltzmann constant, J/K
    T = 300.0; % Temperature, K
    thermal = k_B*T/q; % Thermal voltage, V
    x = Deltax*transpose([0:N-1]); % real space, m
    eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
    Ndon = 2e21*ones(N,1); % 2e15 /cm^3
    Ndon(1:interface1,1) = 5e23; % 5e17 /cm^3
    Ndon(interface2:N,1) = 5e23; % 5e17 /cm^3
    ni = 1.075e16; % 1.075e10 /cm^3
    coef = Deltax*Deltax*q/eps0;
    
    
    
    
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
    
    elec = zeros(N,1);    
    elec = ni*exp(phi/thermal);
    
    if r == 1
        plot(x*1e9,elec/1e6,'r');
        hold on;
    elseif r == 2
        plot(x*1e9,elec/1e6,'b');
        hold on;
    else
        plot(x*1e9,elec/1e6,'g');
        hold on;
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
        res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
        Jaco(2*N-1,2*N-1) = 1.0;
        
        %%% Continuity %%%
        for ii=2:N-1 
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
    
    if r == 1
        plot(x*1e9,elec/1e6,'ro');
        hold on;
    elseif r == 2
        plot(x*1e9,elec/1e6,'bo');
        hold on;
    else
        plot(x*1e9,elec/1e6,'go');
        hold on;
    end
end


xlabel('Position (nm)')
ylabel('electron density (cm^{-3})')

