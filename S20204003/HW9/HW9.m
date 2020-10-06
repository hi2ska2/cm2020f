% Mesh information
Nxx=9;          % X axis mesh
Nyy=5;          % Y axis mesh

%index and voltage application
Source=2; Drain=8; 
Vd=0;           %Drain Bias
Vs=0;            %Source Bias
Vsub=1;         %Substrate Bias
A=sparse(Nxx*Nyy,Nxx*Nyy);

%Matrix
for iy=1:Nyy    %y방향을 행 index로 하고 x방향을 열의 index로 하였습니다(좌표평면과 비슷)
    
    for ix=1+(iy-1)*Nxx:Nxx+(iy-1)*Nxx    %x축에 대한 행렬 구성이 끝나면 그 다음 y index에 대해 계산합니다
    
        if iy==1   %Interace with source and drain
            if ix<=Source     % Dirichlet boundary
            A(ix,ix)=1;  b(ix,1)=Vs;   
            elseif ix>=Drain  % Dirichlet boundary
            A(ix,ix)=1;  b(ix,1)=Vd;      
            else % Neumann boundary
            A(ix,ix)=-2; A(ix,ix+1)=0.5; A(ix,ix-1)=0.5; A(ix,ix+Nxx)=1;
            b(ix,1)=0;
            end
            
        elseif iy==Nyy   %substrate
            A(ix,ix)=1; b(ix,1)=Vsub;
            
        else   %Bulk
            if ix==1+Nxx*(iy-1) % Neumann boundary
            A(ix,ix)=-2; A(ix,ix+1)=1; A(ix,ix+Nxx)=0.5; A(ix,ix-Nxx)=0.5;
            b(ix,1)=0;
            elseif ix==Nxx+Nxx*(iy-1) % Neumann boundary
            A(ix,ix)=-2; A(ix,ix-1)=1; A(ix,ix+Nxx)=0.5; A(ix,ix-Nxx)=0.5;
            b(ix,1)=0;
            else
            A(ix,ix)=-4; A(ix,ix+1)=1; A(ix,ix-1)=1; A(ix,ix+Nxx)=1; A(ix,ix-Nxx)=1;
            b(ix,1)=0;
            end
      
        end
    end
end
% Potential
phi= A\b;
phi=transpose(reshape(phi,[Nxx,Nyy]));

%Plot
x=0:8;
y=4:-1:0;
surf(x,y,phi);
