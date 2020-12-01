clear;

vertex = load('MERGED.vertex.txt');
element = load('MERGED.element.txt');

for row = 1:size(element,1);
    i = [element(row,1:3) element(row,1)];
    p = patch(vertex(i,1),vertex(i,2),'w','LineStyle','-');
end

phi = zeros(size(vertex,1),1);
Jaco = sparse(size(vertex,1),size(vertex,1));
b = zeros(size(vertex,1),1);

for i=1:size(element)
    for ii=0:2         
        vertex_p1 = element(i,rem(ii,3)+1);
        vertex_p2 = element(i,rem(ii+1,3)+1);
        vertex_p3 = element(i,rem(ii+2,3)+1);

        x1 = vertex(vertex_p1,1);
        y1 = vertex(vertex_p1,2);

        x2 = vertex(vertex_p2,1);
        y2 = vertex(vertex_p2,2);

        x3 = vertex(vertex_p3,1);
        y3 = vertex(vertex_p3,2);
        
        p1 = [x1 y1 0]; p2 = [x2 y2 0]; p3 = [x3 y3 0];
        v1 = p1-p2; v2 = p2-p3; v3 = p3-p1;
%%% Length % Area
        L1 = sqrt(v1*v1');
        L2 = sqrt(v2*v2');
        L3 = sqrt(v3*v3');
        
        L = [L1 L2 L3];
        L = sort(L);

        Area = sqrt(cross(v1,v2)*cross(v1,v2)')/2;

        R = L1*L2*L3/4/Area;

        Area12 = sqrt(abs(R^2-(L1/2)^2));
        Area23 = sqrt(abs(R^2-(L2/2)^2));
        Area31 = sqrt(abs(R^2-(L3/2)^2));

%%% Triangle type   
        if (L(1,3)==L2 || L(1,3)^2< L(1,2)^2+L(1,1)^2)
            Area23 = -Area23;
        elseif (L(1,3)==L3 || L(1,3)^2< L(1,2)^2+L(1,1)^2)
            Area31 = -Area31;
        elseif (L(1,3)==L1 || L(1,3)^2< L(1,2)^2+L(1,1)^2)
            Area12 = -Area12;
        end

%%% Jacobian
        Jaco(vertex_p1,vertex_p1) = Jaco(vertex_p1,vertex_p1) - 1/L1*Area12;
        Jaco(vertex_p1,vertex_p2) = Jaco(vertex_p1,vertex_p2) + 1/L1*Area12;
        Jaco(vertex_p1,vertex_p1) = Jaco(vertex_p1,vertex_p1) - 1/L3*Area31;
        Jaco(vertex_p1,vertex_p3) = Jaco(vertex_p1,vertex_p3) + 1/L3*Area31;
    end

end

%%% Boundary condition
Jaco(100,:) = 0; Jaco(100,100) = 1; b(100,1) = 1;
Jaco(500,:) = 0; Jaco(500,500) = 1; b(500,1) = 0;

update = Jaco\b;
phi = phi + update;

figure(2)

for j = 1:size(element,1)
    i = [element(j,1:3)];      
    fill3(vertex(i,1)/1e-9, vertex(i,2)/1e-9,phi(i,1),phi(i,1)); hold on;
end
xlabel('x (nm)');
ylabel('y (nm)');


