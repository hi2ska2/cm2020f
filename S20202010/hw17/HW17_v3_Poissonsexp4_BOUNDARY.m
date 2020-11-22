%Assignment17 initial guess
clear all;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, m
tsi = 5e-9; %silicon layer thickness, m
k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.0e16; % 1.0e10/cm^3, intrinsic carrier density
T=300; %temp. 300K
Nacc = 1e26; %number density, m^-3, 10^20cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
phi0s = 0.5953; % source and drain contact boundary condition
N = 71*301;
phi0g = 0.3374; % gate contact
Jaco1 = sparse(N,N);
phi1 = zeros(N,1);
res1 = zeros(N,1);
nx = 301;
deltax = 1e-10; %0.1nm
Deltax = deltax;
coef = deltax*deltax*q/eps0;
Vg = 0; %gate voltage;
eps_m=(eps_ox+eps_si)/2;
un = 1500e-4;
diajaco = zeros(N,1);

%{
본인 phi((iy-1)*nx+ix,1)     Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)
오른 phi((iy-1)*nx+ix+1,1)   Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)
아래 phi((iy-2)*nx+ix,1)      Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)
왼   phi((iy-1)*nx+ix-1,1)   Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)
위   phi((iy)*nx+ix,1)     Jaco((iy-1)*nx+ix,(iy)*nx+ix)

%}

Ndop = ni*ones(N,1);
% Nonlinear Poisson's

%%%%% doping concentration %%%%%%

for iy = 11:60 %doped layer
    for ix = 1:100
        Ndop((iy-1)*nx+ix,1) = Nacc;
    end
    for ix = 201:301
        Ndop((iy-1)*nx+ix,1) = Nacc;
    end
end



phi = zeros(N,1);
phi(:,1) = V_T*log(Ndop(:,1)/ni);
elec = zeros(N,1);
elec = ni*exp(phi/V_T);
res = zeros(N,1);
Jaco = sparse(N,N);

figure(1);
psi_2d = reshape(phi, 301, 71);
surf(psi_2d);


for nt = 1 :20 %inital non-linear poisson's eq.




        for iy = 1:10 %oxide layer
            for ix = 1:301
                if (iy==1) %아래 첫줄
                    if (ix == 1) %아래 왼쪽 모서리
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 본인
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; %본인 오른 위
                    elseif (ix == nx) %아래 오른쪽 모서리
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %왼 위 본인
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; % 본인 왼 위
                    elseif (ix<101 || ix >201)
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %왼 위 본인 오른
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % 본인 왼 위 오른
                    else %bottom gate 영역
                        res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                    end
                elseif (iy<10)
                        if(ix==1) %왼쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 오른 위 아래 본인
                        elseif (ix == nx) %오른쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 왼 위 아래 본인
                        else
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % 왼 위 아래 오른 본인
                        end
                else %oxide 와 실리콘의경계
                        if(ix==1) %왼쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 오른 위 아래 본인
                        elseif (ix == nx) %오른쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 왼 위 아래 본인
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                        end
                end
            end
        end
        for iy = 11:60 %silicon layer
            for ix = 1:301

                if(iy==11) %oxide와 실리콘의 경계
                    if(ix==1) % 왼쪽 모서리Souce 영역
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %오른쪽 모서리(drain)
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                    end
                elseif(iy<60) %silicon bulk
                    if(ix==1) % 왼쪽 모서리
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %오른쪽 모서리
                          res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                          Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_si*phi((iy-1)*nx+ix-1,1)+eps_si*phi((iy)*nx+ix,1)+eps_si*phi((iy-2)*nx+ix,1)+eps_si*phi((iy-1)*nx+ix+1,1)-4*eps_si*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_si; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_si; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_si; % 왼 위 아래 오른 본인
                    end
                else % intersection between sio2 and si

                    if(ix==1) % 왼쪽 모서리
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %오른쪽 모서리
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                    end
                end
            end
        end



        for iy = 61:71 %top oxide layer
                for ix = 1:nx
                if (iy==71) %맨 윗줄
                    if (ix == 1) %위 왼쪽 모서리
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % 오른아래 본인%%%%%%%%
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; %본인 오른아래
                    elseif (ix == nx) %아래 오른쪽 모서리
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %왼 아래 본인
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; % 본인 왼아래
                    elseif (ix<101 || ix >201)
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %왼 아래 본인 오른
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % 본인 왼아래 오른
                    else %top gate 영역
                        res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                    end %%여기까지함

                elseif (iy == 61) %oxide 와 실리콘의경계
                        if(ix==1) %왼쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 오른 위 아래 본인
                        elseif (ix == nx) %오른쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 왼 위 아래 본인
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                        end
                else % oxide bulk 
                        if(ix==1) %왼쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 오른 위 아래 본인
                        elseif (ix == nx) %오른쪽 모서리
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 왼 위 아래 본인
                        else
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % 왼 위 아래 오른 본인
                        end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%charge part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for iy = 11:60
            for ix = 1: 301
                if(iy==11 || iy == 60) 
                    if (ix==1) % 왼쪽 모서리부분
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.25;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.25;
                    elseif (ix==301)
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.25;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.25;
                    elseif(ix<101 || ix>200)
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                    else
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                    end
                else
                    if (ix==1) % 왼쪽 모서리부분
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                    elseif (ix==301)
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                    elseif (ix<101 || ix>200)
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T));
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T;
                    else
                        res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T));
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T;
                    end
                end
            end
        end 

        update = Jaco\(-res);
        phi(:,1) = phi(:,1)+update ;


end

res4 = zeros(N,1);
res4= res;
    
    figure(2);

    psi_2d=reshape(phi,301,71);
    surf(psi_2d);
        
    Ix = zeros(11,11);
    Vx = zeros(11,1);
    phi1 = phi ;
    %%% scaling %%%
    %%%% Bernouill equation %%%%
    %% S-G scheme
    
    %{
B(x)    
(dphidx/V_T)/(exp(dphidx/V_T)-1)
B(-x)
((-dphidx/V_T)/(exp(-dphidx/V_T) -1))
B'(x)
    (exp(dphidx/V_T)*(dphidx/V_T-1)+1)/(V_T*(exp(dphidx/V_T)-1)^2)
    
    B'(-x)
    (exp((-dphidx/V_T))*((dphidx/V_T)+1)-1)/((-V_T)*(exp(-dphidx/V_T)-1)^2)
    
    
    %}

for bias = 1:11
    for gate = 1:1
       V_applied =  0.05*(bias-1);
       Vx(bias,1) = V_applied;
       Vg = (gate-1)*0.1;
       elec = ni*exp(phi/V_T);
        for nt = 1:5
            res2 = zeros(2*N,1);
            Jaco2 = sparse(2*N,2*N);
            res = zeros(N,1);
            Jaco = sparse(N,N);
            
            
       %%Poisson's eq. %%%
           for iy = 1:10 %oxide layer
                for ix = 1:301
                    if (iy==1) %아래 첫줄
                        if (ix == 1) %아래 왼쪽 모서리
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 본인
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; %본인 오른 위
                        elseif (ix == nx) %아래 오른쪽 모서리
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %왼 위 본인
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; % 본인 왼 위
                        elseif (ix<101 || ix >201)
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %왼 위 본인 오른
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % 본인 왼 위 오른
                        else %bottom gate 영역
                            res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                        end
                    elseif (iy<10)
                            if(ix==1) %왼쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 오른 위 아래 본인
                            elseif (ix == nx) %오른쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 왼 위 아래 본인
                            else
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % 왼 위 아래 오른 본인
                            end
                    else %oxide 와 실리콘의경계
                            if(ix==1) %왼쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 오른 위 아래 본인
                            elseif (ix == nx) %오른쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 왼 위 아래 본인
                            else
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                            end
                    end
                end
            end
            for iy = 11:60 %silicon layer
                for ix = 1:301

                    if(iy==11) %oxide와 실리콘의 경계
                        if(ix==1) % 왼쪽 모서리Souce 영역
                                 res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        elseif (ix==nx) %오른쪽 모서리(drain)
                                 res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                        end
                    elseif(iy<60) %silicon bulk
                        if(ix==1) % 왼쪽 모서리
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        elseif (ix==nx) %오른쪽 모서리
                              res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                              Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        else
                             res((iy-1)*nx+ix,1) = eps_si*phi((iy-1)*nx+ix-1,1)+eps_si*phi((iy)*nx+ix,1)+eps_si*phi((iy-2)*nx+ix,1)+eps_si*phi((iy-1)*nx+ix+1,1)-4*eps_si*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_si; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_si; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_si; % 왼 위 아래 오른 본인
                        end
                    else % intersection between sio2 and si

                        if(ix==1) % 왼쪽 모서리
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        elseif (ix==nx) %오른쪽 모서리
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s-V_applied; % Bias Voltage
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                        end
                    end
                end
            end



            for iy = 61:71 %top oxide layer
                    for ix = 1:nx
                    if (iy==71) %맨 윗줄
                        if (ix == 1) %위 왼쪽 모서리
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % 오른아래 본인%%%%%%%%
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; %본인 오른아래
                        elseif (ix == nx) %아래 오른쪽 모서리
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %왼 아래 본인
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; % 본인 왼아래
                        elseif (ix<101 || ix >201)
                            res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %왼 아래 본인 오른
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % 본인 왼아래 오른
                        else %top gate 영역
                            res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                        end %%여기까지함

                    elseif (iy == 61) %oxide 와 실리콘의경계
                            if(ix==1) %왼쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 오른 위 아래 본인
                            elseif (ix == nx) %오른쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % 왼 위 아래 본인
                            else
                                 res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % 왼 위 아래 오른 본인
                            end
                    else % oxide bulk 
                            if(ix==1) %왼쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 오른 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 오른 위 아래 본인
                            elseif (ix == nx) %오른쪽 모서리
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % 왼 위 아래 본인
                            else
                                 res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % 왼 위 아래 오른 본인
                                 Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % 왼 위 아래 오른 본인
                            end
                    end
                end
            end
       
           
                    

            %%%%%%%%%%%%%%%%%%%%%%%%%%%charge part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for iy = 11:60
                for ix = 1: 301
                    if(iy==11 || iy == 60) 
                        if (ix==1) % 왼쪽 모서리부분
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.25;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.25;
                        elseif (ix==301)
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.25;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.25;
                        elseif(ix<101 || ix>200)
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                        else
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                        end
                    else
                        if (ix==1) % 왼쪽 모서리부분
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                        elseif (ix==301)
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))*0.5;
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T*0.5;
                        elseif (ix<101 || ix>200)
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T));
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(Nacc+ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T;
                        else
                            res((iy-1)*nx+ix,1) = res((iy-1)*nx+ix,1) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T));
                            Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) - coef*(ni*exp(phi((iy-1)*nx+ix,1)/V_T))/V_T;
                        end
                    end
                end
            end 
            for ix = 1:301
                for iy = 1:71
                    xx2 = (iy-1)*nx+ix+1;
                    xx1 = (iy-1)*nx+ix;
                    xx_1 = (iy-1)*nx+ix-1;
                    yy2 = (iy)*nx+ix;
                    yy1 = (iy-1)*nx+ix;
                    yy_1 = (iy-2)*nx+ix;
                    
                    if ix == 301
                        dphidx = 2;
                    else
                        dphidx = phi(xx2,1) - phi(xx1,1);
                    end
                        x1 = dphidx/V_T;
                        cx1 = abs(x1);

                        y=x1;
                        y2 = -x1;
                        if (cx1 < 2.502e-02)

                            sxd = y*y;
                            Bx1 = (1.0 - y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            Bx11 = (1.0 + y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            dBdx1= (-0.5+y/6*(1-sxd/30*(1-sxd/28)));
                            dBdx11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28)));

                        elseif (cx1 < 1.50000e-1)
                            sxd = y*y;
                            Bx1 = (1.0 - y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            Bx11 = (1.0 + y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            dBdx1 = (-0.5+y/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                            dBdx11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                        elseif (cx1 > 150.01)
                            inv_expx = exp(-y);
                            inv_expx2 = exp(-y2);
                            Bx1 = y*inv_expx;
                            Bx11 = y2*inv_expx2;
                            dBdx1 = (inv_expx - Bx1);
                            dBdx11 = (inv_expx2 - Bx11);

                        else
                            inv_expx_1 = 1/(exp(y)-1.0);
                            inv_expx_11 = 1/(exp(y2)-1.0);
                            Bx1 = y*inv_expx_1;
                            Bx11 = y2*inv_expx_11;
                            dBdx1 = (inv_expx_1 - Bx1*(inv_expx_1+1));
                            dBdx11 = (inv_expx_11 - Bx11*(inv_expx_11+1));
                        end
                    if ix == 1
                        dphidx_1 = 2;
                    else
                        dphidx_1 = phi(xx1,1) - phi(xx_1,1);
                    end
                        x_1 = dphidx_1/V_T;
                        cx_1 = abs(x_1);
                        y=x_1;
                        y2 = -x_1;
                          if (cx_1 < 2.502e-02)

                            sxd = y*y;
                            Bx_1 = (1.0 - y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            Bx_11 = (1.0 + y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            dBdx_1= (-0.5+y/6*(1-sxd/30*(1-sxd/28)));
                            dBdx_11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28)));

                        elseif (cx_1 < 1.50000e-1)
                            sxd = y*y;
                            Bx_1 = (1.0 - y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            Bx_11 = (1.0 + y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            dBdx_1 = (-0.5+y/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                            dBdx_11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                        elseif (cx_1 > 150.01)
                            inv_expx = exp(-y);
                            inv_expx2 = exp(-y2);
                            Bx_1 = y*inv_expx;
                            Bx_11 = y2*inv_expx;
                            dBdx_1 = (inv_expx - Bx_1);
                            dBdx_11 = (inv_expx2 - Bx_11);
                        else
                            inv_expx_1 = 1/(exp(y)-1.0);
                            inv_expx_11 = 1/(exp(y2)-1.0);
                            Bx_1 = y*inv_expx_1;
                            Bx_11 = y2*inv_expx_11;
                            dBdx_1 = (inv_expx_1 - Bx_1*(inv_expx_1+1));
                            dBdx_11 = (inv_expx_11 - Bx_11*(inv_expx_11+1));
                          end
                        if iy == 71
                            dphidy = 2;
                        else
                            dphidy = phi(yy2,1) - phi(yy1,1);
                        end
                        y1 = dphidy/V_T;
                        cy1 = abs(y1);
                        y=y1;
                        y2=-y1;
                        if (cy1 < 2.502e-02)

                            sxd = y*y;
                            By1 = (1.0 - y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            By11 = (1.0 + y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            dBdy1= (-0.5+y/6*(1-sxd/30*(1-sxd/28)));
                            dBdy11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28)));

                        elseif (cy1 < 1.50000e-1)
                            sxd = y*y;
                            By1 = (1.0 - y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            By11 = (1.0 + y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            dBdy1 = (-0.5+y/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                            dBdy11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                        elseif (cy1 > 150.01)
                            inv_expx = exp(-y);
                            By1 = y*inv_expx;
                            inv_expx2 = exp(-y2);
                            By11 = y2*inv_expx;
                            dBdy1 = (inv_expx - By1);
                            dBdy11 = (inv_expx2 - By11);
                        else
                            inv_expx_1 = 1/(exp(y)-1.0);
                            By1 = y*inv_expx_1;
                            inv_expx_11 = 1/(exp(y2)-1.0);
                            By11 = y2*inv_expx_11;
                            dBdy1 = (inv_expx_1 - By1*(inv_expx_1+1));
                            dBdy11 = (inv_expx_11 - By11*(inv_expx_11+1));
                        end
                    if iy == 1
                        dphidy_1 = 2 ;
                    else
                        dphidy_1 = phi(yy1,1) - phi(yy_1,1);
                    end
                        y_1 = dphidy/V_T;
                        cy_1 = abs(y_1);
                        y=y_1;
                        y2= -y_1;
                        if (cy_1 < 2.502e-02)

                            sxd = y*y;
                            By_1 = (1.0 - y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            By_11 = (1.0 + y/2 + sxd/12*(1-sxd/60*(1-sxd/42)) );
                            dBdy_1= (-0.5+y/6*(1-sxd/30*(1-sxd/28)));
                            dBdy_11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28)));
                        elseif (cy_1 < 1.50000e-1)
                            sxd = y*y;
                            By_1 = (1.0 - y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            By_11 = (1.0 + y/2.0+sxd/12*(1-sxd/60*(1-sxd/42*(1-sxd/40*(1-0.02525252525252525252525*sxd)))));
                            dBdy_1 = (-0.5+y/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                            dBdy_11 = (-0.5+y2/6*(1-sxd/30*(1-sxd/28*(1-sxd/30*(1-0.031565656565656565656567*sxd)))));
                        elseif (cx1 > 150.01)
                            inv_expx = exp(-y);
                            By_1 = y*inv_expx;
                            inv_expx2 = exp(-y2);
                            By_11 = y2*inv_expx;
                            dBdy_1 = (inv_expx - By_1);
                            dBdy_11 = (inv_expx2 - By_11);
                        else
                            inv_expx_1 = 1/(exp(y)-1.0);
                            By_1 = y*inv_expx_1;
                            inv_expx_11 = 1/(exp(y2)-1.0);
                            By_11 = y2*inv_expx_11;
                            dBdy_1 = (inv_expx_1 - By_1*(inv_expx_1+1));
                            dBdy_11 = (inv_expx_11 - By_11*(inv_expx_11+1));
                        end
                    
                    
                    
                    
                    if (iy==1)
                        if (ix == 1) %왼쪽 아래 모서리 오른 위만 존재
                            res2(2*xx1,1) = 0.5*(elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11) + 0.5*(elec(yy2,1)*By1-elec(yy1,1)*By11);
                            Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + 0.5*(elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);
                            Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + 0.5*(elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T);

                            Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + 0.5*(elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T + elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T);
                            Jaco2(2*xx1, 2*xx2) =Jaco2(2*xx1, 2*xx2)+ 0.5*Bx1;
                            Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2)+0.5*By1;

                            Jaco2(2*xx1, 2*xx1) =Jaco2(2*xx1, 2*xx1)+ 0.5*(-Bx11 - By11);    
                        elseif (ix == 301) % 오른아래 모서리 (왼 위 만)
                            res2(2*xx1,1) = 0.5*(elec(yy2,1)*By1-elec(yy1,1)*By11)  -   (  0.5*(elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11) ); %오른 위  왼 아래
                          
                            Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + 0.5*(elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T); %위
                            Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + 0.5*(elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                            
                            Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + 0.5*(elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T) +0.5*(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) ;
                            
                            Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2) + 0.5*By1; % 위
                            Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) +0.5*Bx_11; %왼
                            
                            Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1)  - 0.5*By11 - 0.5*Bx_1 ; %오른 위 왼 아래

                        else %아래 한줄 (0.5오른 위 0.5왼 )
                            res2(2*xx1,1) = 0.5*(elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11) + (elec(yy2,1)*By1-elec(yy1,1)*By11)  -   0.5*(  (elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11)) ; %오른 위  왼 아래
                            Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + 0.5*(elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);%오른
                            Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + (elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T); %위
                            Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + 0.5*(elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                            
                            Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + 0.5*(elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T) + (elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T) +0.5*(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) ;
                            Jaco2(2*xx1, 2*xx2) = Jaco2(2*xx1, 2*xx2) + 0.5*Bx1; % 오른
                            Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2) + By1; % 위
                            Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) + 0.5*Bx_11; %왼
                           
                            Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1) -0.5*Bx11 - By11 - 0.5*Bx_1; %오른 위 왼 아래
                  
                            
                            
                        end
                    elseif (iy == 71)
                        
                            if (ix == 1) %왼쪽 위 모서리 (오른 아래만)
                                res2(2*xx1,1) = 0.5*(elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11  - (elec(yy1,1)*By_1-elec(yy_1,1)*By_11));
                                Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + 0.5*(elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);%오른
                                
                               
                                Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + 0.5*(elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                                Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + 0.5*(elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T - elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                                Jaco2(2*xx1, 2*xx2) = 0.5*Bx1; % 오른
                               
                                Jaco2(2*yy1, 2*yy_1) = 0.5*By_11; %아래
                                Jaco2(2*xx1, 2*xx1) = -0.5*Bx11 - 0.5*By_1; %오른 위 왼 아래
                            
                            elseif (ix == 301) % 오른 위 모서리 (왼아래만)
                                res2(2*xx1,1) =  -   0.5*(  (elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11)  +   (elec(yy1,1)*By_1-elec(yy_1,1)*By_11) ); %오른 위  왼 아래
                                
                                Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + 0.5*(elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                                Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + 0.5*(elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                                Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1)  +0.5*(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) +0.5*(- elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                              
                                Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) + 0.5*Bx_11; %왼
                                Jaco2(2*yy1, 2*yy_1) = Jaco2(2*yy1, 2*yy_1) + 0.5*By_11; %아래
                                Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1)  - 0.5*Bx_1 - 0.5*By_1; %오른 위 왼 아래
                            else %위 한줄 (아래1 왼 오른 반)
                                res2(2*xx1,1) = 0.5*(elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11)   -   (  0.5*(elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11)  +   (elec(yy1,1)*By_1-elec(yy_1,1)*By_11) ); %오른   왼 아래
                                Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + 0.5*(elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);%오른
                                 
                                Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + 0.5*(elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                                Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + (elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                                Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + 0.5*(elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T)  +0.5*(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) +(- elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                                Jaco2(2*xx1, 2*xx2) = Jaco2(2*xx1, 2*xx2) + 0.5*Bx1; % 오른
                            
                                Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) + 0.5*Bx_11; %왼
                                Jaco2(2*yy1, 2*yy_1) = Jaco2(2*yy1, 2*yy_1) + By_11; %아래
                                Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1) -0.5*Bx11 - 0.5*Bx_1 - By_1; %오른 위 왼 아래
                                
                            end
                            
                        
                    else %bulk
                        if (ix ==1) %왼 한줄 (왼 없음) 위 아래 반
                            res2(2*xx1,1) = (elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11) + 0.5*(elec(yy2,1)*By1-elec(yy1,1)*By11)  -   0.5*((elec(yy1,1)*By_1-elec(yy_1,1)*By_11) ); %오른 위  왼 아래
                            Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + (elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);%오른
                            Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + 0.5*(elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T); %위
                         
                            Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + 0.5*(elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                            Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + (elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T) + 0.5*(elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T)  +0.5*(- elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                            Jaco2(2*xx1, 2*xx2) = Jaco2(2*xx1, 2*xx2) + Bx1; % 오른
                            Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2) + 0.5*By1; % 위
                            
                            Jaco2(2*yy1, 2*yy_1) = Jaco2(2*yy1, 2*yy_1) + 0.5*By_11; %아래
                            Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1) -Bx11 - 0.5*By11  - 0.5*By_1; %오른 위 왼 아래
                        elseif (ix == 301) %오른한줄
                            res2(2*xx1,1) = (elec(yy2,1)*By1-elec(yy1,1)*By11)  -   (  0.5*(elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11)  +   0.5*(elec(yy1,1)*By_1-elec(yy_1,1)*By_11) ); %오른 위  왼 아래
                           
                            Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + 0.5*(elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T); %위
                            Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + (elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                            Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + 0.5*(elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                            Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) +  0.5*(elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T) +(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) +0.5*(- elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                            
                            Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2) + 0.5*By1; % 위
                            Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) + Bx_11; %왼
                            Jaco2(2*yy1, 2*yy_1) = Jaco2(2*yy1, 2*yy_1) + 0.5*By_11; %아래
                            Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1) - 0.5*By11 - Bx_1 - 0.5*By_1; %오른 위 왼 아래
                        
                        else
                        res2(2*xx1,1) = (elec(xx2,1)*Bx1 - elec(xx1,1)*Bx11) + (elec(yy2,1)*By1-elec(yy1,1)*By11)  -   (  (elec(xx1,1)*Bx_1-elec(xx_1,1)*Bx_11)  +   (elec(yy1,1)*By_1-elec(yy_1,1)*By_11) ); %오른 위  왼 아래
                        Jaco2(2*xx1, 2*xx2-1) = Jaco2(2*xx1, 2*xx2-1) + (elec(xx2,1)*dBdx1/V_T + elec(xx1,1)*dBdx11/V_T);%오른
                        Jaco2(2*yy1, 2*yy2-1) = Jaco2(2*yy1, 2*yy2-1) + (elec(yy2,1)*dBdy1/V_T + elec(yy1,1)*dBdy11/V_T); %위
                        Jaco2(2*xx1, 2*xx_1-1) = Jaco2(2*xx1, 2*xx_1-1) + (elec(xx1,1)*dBdx_1/V_T + elec(xx_1,1)*dBdx_11/V_T); % 왼
                        Jaco2(2*yy1, 2*yy_1-1) = Jaco2(2*yy1, 2*yy_1-1) + (elec(yy1,1)*dBdy_1/V_T + elec(yy_1,1)*dBdy_11/V_T); % 아래
                        Jaco2(2*xx1, 2*xx1-1) = Jaco2(2*xx1, 2*xx1-1) + (elec(xx2,1)*dBdx1/(-V_T) - elec(xx1,1)*dBdx11/V_T) + (elec(yy2,1)*dBdy1/(-V_T) - elec(yy1,1)*dBdy11/V_T) +(- elec(xx1,1)*dBdx_1/V_T+elec(xx_1,1)*dBdx_11/(-V_T)) +(- elec(yy1,1)*dBdy_1/V_T+elec(yy_1,1)*dBdy_11/(-V_T));
                        Jaco2(2*xx1, 2*xx2) = Jaco2(2*xx1, 2*xx2) + Bx1; % 오른
                        Jaco2(2*yy1, 2*yy2) = Jaco2(2*yy1, 2*yy2) + By1; % 위
                        Jaco2(2*xx1, 2*xx_1) = Jaco2(2*xx1, 2*xx_1) + Bx_11; %왼
                        Jaco2(2*yy1, 2*yy_1) = Jaco2(2*yy1, 2*yy_1) + By_11; %아래
                        Jaco2(2*xx1, 2*xx1) = Jaco2(2*xx1, 2*xx1) -Bx11 - By11 - Bx_1 - By_1; %오른 위 왼 아래
                        
                        end 
                    end
                end
            end
      
                    
                    
                    
                    
            
                     
            
            for ix = 1:301
                for iy = 1:71

                    res2(2*((iy-1)*nx+ix)-1,1) = res(((iy-1)*nx+ix),1);

                    Jaco2(2*((iy-1)*nx+ix)-1,2*((iy-1)*nx+ix)-1) = Jaco(((iy-1)*nx+ix),((iy-1)*nx+ix));
                    if(ix<301)
                        Jaco2(2*((iy-1)*nx+ix)-1,2*((iy-1)*nx+ix+1)-1) = Jaco(((iy-1)*nx+ix),((iy-1)*nx+ix+1));
                    end
                    if(iy<71)
                        Jaco2(2*((iy-1)*nx+ix)-1,2*((iy)*nx+ix)-1) = Jaco(((iy-1)*nx+ix),((iy)*nx+ix));
                    end
                    if(ix>1)
                        Jaco2(2*((iy-1)*nx+ix)-1,2*((iy-1)*nx+ix-1)-1) = Jaco(((iy-1)*nx+ix),((iy-1)*nx+ix-1));
                    end
                    if(iy>1)
                        Jaco2(2*((iy-1)*nx+ix)-1,2*((iy-2)*nx+ix)-1) = Jaco(((iy-1)*nx+ix),((iy-2)*nx+ix));
                    end
                end
            end
               
            for ix = 2:300
                for iy = 2:70
                    Jaco2(2*((iy-1)*nx+ix)-1,2*((iy-1)*nx+ix)) = -coef;
                end
            end
            
          for ix = 1:301
            for iy = 1:71
                if(iy == 1 || iy ==71)
                    if(ix >101 || ix<201)
                        res2(2*((iy-1)*nx+ix), 1) = elec(((iy-1)*nx+ix), 1) - Ndop(((iy-1)*nx+ix),1);
                        Jaco2(2*((iy-1)*nx+ix), :) = 0.0;
                        Jaco2(2*((iy-1)*nx+ix), 2*((iy-1)*nx+ix)) = 1.0;
                    end
                elseif(iy >10 || iy<61)
                    if(ix == 1 ||ix==301)
                        res2(2*((iy-1)*nx+ix), 1) = elec(((iy-1)*nx+ix), 1) - Ndop(((iy-1)*nx+ix),1);
                        Jaco2(2*((iy-1)*nx+ix), :) = 0.0;
                        Jaco2(2*((iy-1)*nx+ix), 2*((iy-1)*nx+ix)) = 1.0;

                    end
                end
            end
        end
            
       if (nt ==1)
           for i = 1:N
               diajaco(i,1) = Jaco2(i,i);
           end
       end
            
        Cvector = zeros(2*N,1);
        Cvector(1:2:2*N-1,1) = V_T;
        Cvector(2:2:2*N,1) = Nacc;
        Cmatrix = spdiags(Cvector, 0 , 2*N, 2*N);

        Jaco_scaled = Jaco2*Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,2*N,2*N);
        Jaco_scaled = Rmatrix*Jaco_scaled;
        res_scaled = Rmatrix*res2;
        update_scaled = Jaco_scaled\(-res_scaled);
        update = Cmatrix * update_scaled;

        phi = phi + update(1:2:2*N-1,1);
        elec = elec + update(2:2:2*N,1);
            
        end
       for iy = 11:60
           ix = 301;
           Ix(bias, gate) = Ix(bias,gate) + q*un*(elec((iy-1)*nx+ix,1)*(phi((iy-1)*nx+ix,1)-phi((iy-1)*nx+ix-1,1))/Deltax-V_T*(elec((iy-1)*nx+ix,1)-elec((iy-1)*nx+ix-1,1))/Deltax);
       end
    end
end


figure(3)
psi_2d = reshape(Ix/1e9, 11,11);
surf(psi_2d);

      
       
     
           
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 