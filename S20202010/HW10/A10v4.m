%Assignment 10
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
eps_si = 11.7; % relative permittivity for silicon
eps_ox = 3.9;  % relative permittivity for oxide
tox = 8e-10; %oxidelayer thickness, m
tsi = 5e-9; %silicon layer thickness, m
N = 67; % 6.6 nm thick
k_B = 1.380662e-23; % Boltzmann constant, J/K
ni = 1.0e16; % 1.0e10/cm^3, intrinsic carrier density
T=300; %temp. 300K
Nacc = 1e26; %number density, m^-3, 10^20cm^-3
V_T=k_B*T/q; %thermal voltage at 300K (~26meV)
phi0s = 0.5953; % source and drain contact boundary condition
N = 71*301;
phi0g = 0.3374; % gate contact
Jaco = sparse(N,N);
phi = zeros(N,1);
res = zeros(N,1);
nx = 301;
deltax = 1e-10; %0.1nm
coef = deltax*deltax*q/eps0;
Vg = 0; %gate voltage;
eps_m=(eps_ox+eps_si)/2;
%%%%%%%%%%%%%%%%%%%%%%%%% Boundary condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
���� phi((iy-1)*nx+ix,1)     Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)
���� phi((iy-1)*nx+ix+1,1)   Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)
�Ʒ� phi((iy-2)*nx+ix,1)      Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)
��   phi((iy-1)*nx+ix-1,1)   Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)
��   phi((iy)*nx+ix,1)     Jaco((iy-1)*nx+ix,(iy)*nx+ix)

%}

 % initial guess

phi(:,1) = -0.2;

%%%%%%%%%%%%%%% Laplacian part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for z = 1:2
    if z == 2  
        Vg = Vg+1.1;
    end
    for nt = 1 :100




        for iy = 1:10 %oxide layer
            for ix = 1:301
                if (iy==1) %�Ʒ� ù��
                    if (ix == 1) %�Ʒ� ���� �𼭸�
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % ���� �� ����
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; %���� ���� ��
                    elseif (ix == nx) %�Ʒ� ������ �𼭸�
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %�� �� ����
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; % ���� �� ��
                    elseif (ix<101 || ix >201)
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %�� �� ���� ����
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % ���� �� �� ����
                    else %bottom gate ����
                        res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                    end
                elseif (iy<10)
                        if(ix==1) %���� �𼭸�
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % ���� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % ���� �� �Ʒ� ����
                        elseif (ix == nx) %������ �𼭸�
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % �� �� �Ʒ� ����
                        else
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % �� �� �Ʒ� ���� ����
                        end
                else %oxide �� �Ǹ����ǰ��
                        if(ix==1) %���� �𼭸�
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % ���� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % ���� �� �Ʒ� ����
                        elseif (ix == nx) %������ �𼭸�
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % �� �� �Ʒ� ����
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % �� �� �Ʒ� ���� ����
                        end
                end
            end
        end
        for iy = 11:60 %silicon layer
            for ix = 1:301

                if(iy==11) %oxide�� �Ǹ����� ���
                    if(ix==1) % ���� �𼭸�Souce ����
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %������ �𼭸�(drain)
                             res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % �� �� �Ʒ� ���� ����
                    end
                elseif(iy<60) %silicon bulk
                    if(ix==1) % ���� �𼭸�
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %������ �𼭸�
                          res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                          Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_si*phi((iy-1)*nx+ix-1,1)+eps_si*phi((iy)*nx+ix,1)+eps_si*phi((iy-2)*nx+ix,1)+eps_si*phi((iy-1)*nx+ix+1,1)-4*eps_si*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_si; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_si; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_si;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_si; % �� �� �Ʒ� ���� ����
                    end
                else % intersection between sio2 and si

                    if(ix==1) % ���� �𼭸�
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    elseif (ix==nx) %������ �𼭸�
                         res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0s;
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;
                    else
                         res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                         Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % �� �� �Ʒ� ���� ����
                    end
                end
            end
        end



        for iy = 61:71 %top oxide layer
                for ix = 1:nx
                if (iy==71) %�� ����
                    if (ix == 1) %�� ���� �𼭸�
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix+1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); % �����Ʒ� ����%%%%%%%%
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; %���� �����Ʒ�
                    elseif (ix == nx) %�Ʒ� ������ �𼭸�
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*0.5*phi((iy-2)*nx+ix,1)-eps_ox*phi((iy-1)*nx+ix,1); %�� �Ʒ� ����
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; % ���� �޾Ʒ�
                    elseif (ix<101 || ix >201)
                        res((iy-1)*nx+ix,1) = eps_ox*0.5*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1)+0.5*eps_ox*phi((iy-1)*nx+ix+1,1); %�� �Ʒ� ���� ����
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)= -2*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= 0.5*eps_ox; % ���� �޾Ʒ� ����
                    else %top gate ����
                        res((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1)-phi0g-Vg;
                        Jaco((iy-1)*nx+ix,(iy-1)*nx+ix) = 1;

                    end %%���������

                elseif (iy == 61) %oxide �� �Ǹ����ǰ��
                        if(ix==1) %���� �𼭸�
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix+1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % ���� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % ���� �� �Ʒ� ����
                        elseif (ix == nx) %������ �𼭸�
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+0.5*eps_m*phi((iy)*nx+ix,1)+0.5*eps_m*phi((iy-2)*nx+ix,1)-2*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_m; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_m; % �� �� �Ʒ� ����
                        else
                             res((iy-1)*nx+ix,1) = eps_m*phi((iy-1)*nx+ix-1,1)+eps_m*phi((iy)*nx+ix,1)+eps_m*phi((iy-2)*nx+ix,1)+eps_m*phi((iy-1)*nx+ix+1,1)-4*eps_m*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_m; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_m; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_m;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_m; % �� �� �Ʒ� ���� ����
                        end
                else % oxide bulk 
                        if(ix==1) %���� �𼭸�
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix+1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % ���� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % ���� �� �Ʒ� ����
                        elseif (ix == nx) %������ �𼭸�
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+0.5*eps_ox*phi((iy)*nx+ix,1)+0.5*eps_ox*phi((iy-2)*nx+ix,1)-2*eps_ox*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= 0.5*eps_ox; Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-2*eps_ox; % �� �� �Ʒ� ����
                        else
                             res((iy-1)*nx+ix,1) = eps_ox*phi((iy-1)*nx+ix-1,1)+eps_ox*phi((iy)*nx+ix,1)+eps_ox*phi((iy-2)*nx+ix,1)+eps_ox*phi((iy-1)*nx+ix+1,1)-4*eps_ox*phi((iy-1)*nx+ix,1); % �� �� �Ʒ� ���� ����
                             Jaco((iy-1)*nx+ix,(iy-1)*nx+ix-1)= eps_ox; Jaco((iy-1)*nx+ix,(iy)*nx+ix)= eps_ox; Jaco((iy-1)*nx+ix,(iy-2)*nx+ix)= eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix+1)=eps_ox;Jaco((iy-1)*nx+ix,(iy-1)*nx+ix)=-4*eps_ox; % �� �� �Ʒ� ���� ����
                        end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%charge part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for iy = 11:60
            for ix = 1: 301
                if(iy==11 || iy == 60) 
                    if (ix==1) % ���� �𼭸��κ�
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
                    if (ix==1) % ���� �𼭸��κ�
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
    
    figure(z);

    psi_2d=reshape(phi,301,71);
    surf(psi_2d);
end
