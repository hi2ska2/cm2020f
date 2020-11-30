vert=load('MERGED.vertex');
element=load('MERGED.element');

A=sparse(size(vert,1),size(vert,1));
b=zeros(size(vert,1),1);

for ii=1:size(element)
    
    for node=0:2   % To consider all vertexes in an element
%% For element and vertex search
        
        vertex_index1=element(ii,rem(node,3)+1);
        vertex_index2=element(ii,rem(node+1,3)+1);
        vertex_index3=element(ii,rem(node+2,3)+1);
        
        x1 = vert(vertex_index1,1);
        y1 = vert(vertex_index1,2);
        
        x2 = vert(vertex_index2,1);
        y2 = vert(vertex_index2,2);
        
        x3 = vert(vertex_index3,1);
        y3 = vert(vertex_index3,2);
%% For area and length calculation
        % Number 1 is main node
        
        L1 = sqrt((x2-x1)^2 + (y2-y1)^2);
        L2 = sqrt((x3-x2)^2 + (y3-y2)^2);
        L3 = sqrt((x3-x1)^2 + (y3-y1)^2);
        
        Area = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2;
        
        radius=L1*L2*L3/4/Area;
        
        Area12=real(sqrt(radius^2-(L1/2)^2));
        Area23=real(sqrt(radius^2-(L2/2)^2));
        Area31=real(sqrt(radius^2-(L3/2)^2));
        
%% For check triangle type
        
        inner_product1=(x2-x1)*(x3-x1)+(y2-y1)*(y3-y1);
        inner_product2=(x3-x2)*(x1-x2)+(y3-y2)*(y1-y2);
        inner_product3=(x2-x3)*(x1-x3)+(y2-y3)*(y1-y3);
        
        %        patch([x1 x2 x3],[y1 y2 y3],'w','LineStyle','-'); hold on;
        
        if (inner_product1 < 0 || inner_product2 < 0 || inner_product3 < 0 )
            disp(sprintf('%d is obtuse triangle', ii));
            %                 patch([x1 x2 x3],[y1 y2 y3],'r');
        end
        
        if (inner_product1 < 0)
            Area23 = -Area23;
        elseif (inner_product2 < 0)
            Area31 = -Area31;
        elseif (inner_product3 < 0)
            Area12 = -Area12;
        end
        
%% Maxtrix implementation
        A(vertex_index1,vertex_index1) = A(vertex_index1,vertex_index1) - 1/L1*Area12;
        A(vertex_index1,vertex_index2) = A(vertex_index1,vertex_index2) + 1/L1*Area12;
        A(vertex_index1,vertex_index1) = A(vertex_index1,vertex_index1) - 1/L3*Area31;
        A(vertex_index1,vertex_index3) = A(vertex_index1,vertex_index3) + 1/L3*Area31;
    end
    
end
%%%%%%%%Consider Boundary condition%%%%%%%%%%%%%%%%
A(10,:)=0;
A(10,10)=1;
b(10,1)=1;

A(400,:)=0;
A(400,400)=1;
b(400,1)=1;

A(800,:)=0;
A(800,800)=1;
b(800,1)=0;
%%%%%%%%Check Null Matrix Row%%%%%%%%%%%%%%%%%%%%%%
for iii=1:size(vert,1)
    if A (iii,:) == 0
        disp(sprintf('%d Has Null ROW', iii));
    end
end

phi=A \ b;
%% For Plot
vertex_index1=element(:,1); vertex_index2=element(:,2); vertex_index3=element(:,3);
x1 = vert(vertex_index1,1); y1 = vert(vertex_index1,2);
x2 = vert(vertex_index2,1); y2 = vert(vertex_index2,2);
x3 = vert(vertex_index3,1); y3 = vert(vertex_index3,2);
phi1 = phi(vertex_index1,1); phi2 = phi(vertex_index2,1); phi3 = phi(vertex_index3,1);

%  %Type1 : Graph with without edge interpolation
%  patch(transpose([x1 x2 x3]),transpose([y1 y2 y3]),transpose([phi1 phi2 phi3]),'FaceColor','interp');

% % Type2 : Graph with with edge interpolation
%  patch(transpose([x1 x2 x3]),transpose([y1 y2 y3]),transpose([phi1 phi2 phi3]),'FaceColor','interp','EdgeColor','interp');

% %Type3 : 3D Graph without color
% patch(transpose([x1 x2 x3]),transpose([y1 y2 y3]),transpose([phi1 phi2 phi3]),'w','LineStyle','-');
% view(3)

% %Type4 : 3D Graph with color without edge interpolation
% fill3(transpose([x1 x2 x3]),transpose([y1 y2 y3]),transpose([phi1 phi2 phi3]),transpose([phi1 phi2 phi3]));
% view(3)

% Type5 : 3D Graph with color with edge interpolation
fill3(transpose([x1 x2 x3]),transpose([y1 y2 y3]),transpose([phi1 phi2 phi3]),transpose([phi1 phi2 phi3]),'EdgeColor','interp');
view(3)

colorbar
axis square