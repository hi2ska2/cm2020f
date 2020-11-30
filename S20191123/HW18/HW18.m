clear;

element = load('HW18.element');
vertex = load('HW18.vertex');

% for row=1:size(element,1)
%     ii = [element(row,1:3) element(row,1)];
%     p = patch(vertex(ii,1), vertex(ii,2), 'w', 'LineStyle','-');
% end

phi = ones(size(vertex,1),1); % �� vertex������ potential

n1 = 10; %1
n2 = 12; %0


res = zeros(size(vertex,1),1);
Jaco = sparse(size(vertex,1),size(vertex,1));


for newton = 1:5
    %     newton
    for ii = 1:size(element,1)
        element_phi1 = phi(element(ii,1),1); %element�� 1�� vertex�� potential
        element_phi2 = phi(element(ii,2),1); %element�� 2�� vertex�� potential
        element_phi3 = phi(element(ii,3),1); %element�� 3�� vertex�� potential
        
        element_vertex1 = vertex(element(ii,1),:); %element�� 1�� vertex�� ��ǥ
        element_vertex2 = vertex(element(ii,2),:); %element�� 2�� vertex�� ��ǥ
        element_vertex3 = vertex(element(ii,3),:); %element�� 3�� vertex�� ��ǥ
        
        L12 = sqrt((element_vertex1(1,1)-element_vertex2(1,1))^2 + (element_vertex1(1,2)-element_vertex2(1,2))^2 ); %1��2 ���� �Ÿ�
        L23 = sqrt((element_vertex2(1,1)-element_vertex3(1,1))^2 + (element_vertex2(1,2)-element_vertex3(1,2))^2 ); %2��3 ���� �Ÿ�
        L31 = sqrt((element_vertex3(1,1)-element_vertex1(1,1))^2 + (element_vertex3(1,2)-element_vertex1(1,2))^2 ); %3��1 ���� �Ÿ�
        
        L(1,1) = L12;
        L(1,2) = L23;
        L(1,3) = L31;
        
        %����
        Area = 0.5*abs(element_vertex1(1,1)*element_vertex3(1,2)-element_vertex2(1,1)*element_vertex3(1,2)+element_vertex2(1,1)*element_vertex1(1,2)-element_vertex3(1,1)*element_vertex1(1,2)+element_vertex3(1,1)*element_vertex2(1,2)-element_vertex1(1,1)*element_vertex2(1,2));
        
        % �ݰ�
        R = L12*L23*L31/4/Area;
        
        % control area
        A1 = real(sqrt(R^2-(L23/2)^2));
        A2 = real(sqrt(R^2-(L31/2)^2));
        A3 = real(sqrt(R^2-(L12/2)^2));
        
        
        %�а��ﰢ�� �Ǻ�
        if  ((element_vertex1(1,1)-element_vertex2(1,1))*(element_vertex3(1,1)-element_vertex2(1,1)) + (element_vertex1(1,2)-element_vertex2(1,2))*(element_vertex3(1,2)-element_vertex2(1,2)) < 0) || ...
                ((element_vertex2(1,1)-element_vertex3(1,1))*(element_vertex1(1,1)-element_vertex3(1,1)) + (element_vertex2(1,2)-element_vertex3(1,2))*(element_vertex1(1,2)-element_vertex3(1,2)) < 0) ||...
                ((element_vertex3(1,1)-element_vertex1(1,1))*(element_vertex2(1,1)-element_vertex1(1,1)) + (element_vertex3(1,2)-element_vertex1(1,2))*(element_vertex2(1,2)-element_vertex1(1,2)) < 0)
            
            % �а��ﰢ�� �϶�
            if max(L) == L(1,3) %L31
                %             A2(-)
                res(element(ii,1),1) = res(element(ii,1),1) + (-(element_phi3-element_phi1)*A2/L31) + ((element_phi2-element_phi1)*A3/L12);
                res(element(ii,2),1) = res(element(ii,2),1) + ((element_phi3-element_phi2)*A1/L23) + ((element_phi1-element_phi2)*A3/L12);
                res(element(ii,3),1) = res(element(ii,3),1) + (-(element_phi1-element_phi3)*A2/L31) + ((element_phi2-element_phi3)*A1/L23);
                
                Jaco(element(ii,1),element(ii,1)) = Jaco(element(ii,1),element(ii,1)) + (A2/L31) + ((-1)*A3/L12);
                Jaco(element(ii,1),element(ii,2)) = Jaco(element(ii,1),element(ii,2)) + (A3/L12);
                Jaco(element(ii,1),element(ii,3)) = Jaco(element(ii,1),element(ii,3)) + (-A2/L31);
                
                Jaco(element(ii,2),element(ii,1)) = Jaco(element(ii,2),element(ii,1)) + (A3/L12);
                Jaco(element(ii,2),element(ii,2)) = Jaco(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + ((-1)*A3/L12);
                Jaco(element(ii,2),element(ii,3)) = Jaco(element(ii,2),element(ii,3)) + (A1/L23);
                
                Jaco(element(ii,3),element(ii,1)) = Jaco(element(ii,3),element(ii,1)) + (-A2/L31);
                Jaco(element(ii,3),element(ii,2)) = Jaco(element(ii,3),element(ii,2)) + (A1/L23);
                Jaco(element(ii,3),element(ii,3)) = Jaco(element(ii,3),element(ii,3)) + (A2/L31) + ((-1)*A1/L23);
                
            elseif max(L) == L(1,1) %L12
                %             A3(-)
                res(element(ii,1),1) = res(element(ii,1),1) + ((element_phi3-element_phi1)*A2/L31) + (-(element_phi2-element_phi1)*A3/L12);
                res(element(ii,2),1) = res(element(ii,2),1) + ((element_phi3-element_phi2)*A1/L23) + (-(element_phi1-element_phi2)*A3/L12);
                res(element(ii,3),1) = res(element(ii,3),1) + ((element_phi1-element_phi3)*A2/L31) + ((element_phi2-element_phi3)*A1/L23);
                
                Jaco(element(ii,1),element(ii,1)) = Jaco(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + (A3/L12);
                Jaco(element(ii,1),element(ii,2)) = Jaco(element(ii,1),element(ii,2)) + (-A3/L12);
                Jaco(element(ii,1),element(ii,3)) = Jaco(element(ii,1),element(ii,3)) + (A2/L31);
                
                Jaco(element(ii,2),element(ii,1)) = Jaco(element(ii,2),element(ii,1)) + (-A3/L12);
                Jaco(element(ii,2),element(ii,2)) = Jaco(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + (A3/L12);
                Jaco(element(ii,2),element(ii,3)) = Jaco(element(ii,2),element(ii,3)) + (A1/L23);
                
                Jaco(element(ii,3),element(ii,1)) = Jaco(element(ii,3),element(ii,1)) + (A2/L31);
                Jaco(element(ii,3),element(ii,2)) = Jaco(element(ii,3),element(ii,2)) + (A1/L23);
                Jaco(element(ii,3),element(ii,3)) = Jaco(element(ii,3),element(ii,3)) + ((-1)*A2/L31) + ((-1)*A1/L23);
                
            elseif max(L) == L(1,2) %L23
                %             A1(-)
                res(element(ii,1),1) = res(element(ii,1),1) + ((element_phi3-element_phi1)*A2/L31) + ((element_phi2-element_phi1)*A3/L12);
                res(element(ii,2),1) = res(element(ii,2),1) + (-(element_phi3-element_phi2)*A1/L23) + ((element_phi1-element_phi2)*A3/L12);
                res(element(ii,3),1) = res(element(ii,3),1) + ((element_phi1-element_phi3)*A2/L31) + (-(element_phi2-element_phi3)*A1/L23);
                
                Jaco(element(ii,1),element(ii,1)) = Jaco(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + ((-1)*A3/L12);
                Jaco(element(ii,1),element(ii,2)) = Jaco(element(ii,1),element(ii,2)) + (A3/L12);
                Jaco(element(ii,1),element(ii,3)) = Jaco(element(ii,1),element(ii,3)) + (A2/L31);
                
                Jaco(element(ii,2),element(ii,1)) = Jaco(element(ii,2),element(ii,1)) + (A3/L12);
                Jaco(element(ii,2),element(ii,2)) = Jaco(element(ii,2),element(ii,2)) + (A1/L23) + (-A3/L12);
                Jaco(element(ii,2),element(ii,3)) = Jaco(element(ii,2),element(ii,3)) + (-A1/L23);
                
                Jaco(element(ii,3),element(ii,1)) = Jaco(element(ii,3),element(ii,1)) + (A2/L31);
                Jaco(element(ii,3),element(ii,2)) = Jaco(element(ii,3),element(ii,2)) + (-A1/L23);
                Jaco(element(ii,3),element(ii,3)) = Jaco(element(ii,3),element(ii,3)) + (-A2/L31) + (A1/L23);
                
            end
            
        else
            %���� �Ǵ� ���� �ﰢ���� ��
            res(element(ii,1),1) = res(element(ii,1),1) + ((element_phi3-element_phi1)*A2/L31) + ((element_phi2-element_phi1)*A3/L12);
            res(element(ii,2),1) = res(element(ii,2),1) + ((element_phi3-element_phi2)*A1/L23) + ((element_phi1-element_phi2)*A3/L12);
            res(element(ii,3),1) = res(element(ii,3),1) + ((element_phi1-element_phi3)*A2/L31) + ((element_phi2-element_phi3)*A1/L23);
            
            
            Jaco(element(ii,1),element(ii,1)) = Jaco(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + ((-1)*A3/L12);
            Jaco(element(ii,1),element(ii,2)) = Jaco(element(ii,1),element(ii,2)) + (A3/L12);
            Jaco(element(ii,1),element(ii,3)) = Jaco(element(ii,1),element(ii,3)) + (A2/L31);
            
            Jaco(element(ii,2),element(ii,1)) = Jaco(element(ii,2),element(ii,1)) + (A3/L12);
            Jaco(element(ii,2),element(ii,2)) = Jaco(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + ((-1)*A3/L12);
            Jaco(element(ii,2),element(ii,3)) = Jaco(element(ii,2),element(ii,3)) + (A1/L23);
            
            Jaco(element(ii,3),element(ii,1)) = Jaco(element(ii,3),element(ii,1)) + (A2/L31);
            Jaco(element(ii,3),element(ii,2)) = Jaco(element(ii,3),element(ii,2)) + (A1/L23);
            Jaco(element(ii,3),element(ii,3)) = Jaco(element(ii,3),element(ii,3)) + ((-1)*A2/L31) + ((-1)*A1/L23);
        end
        
    end
    
    res(n1,1) = phi(n1,1) - 1;
    Jaco(n1,:) =0;
    Jaco(n1,n1) =1;
    
    res(n2,1) = phi(n2,1);
    Jaco(n2,:) =0;
    Jaco(n2,n2) =1;
    
    update = Jaco\(-res);
    phi = phi +update;
    
end

for row=1:size(element,1)
    ii = [element(row,1:3) element(row,1)];
    
    %     p1 = patch(vertex(ii,1), vertex(ii,2), phi(ii,1),'w', 'LineStyle','-');    
    %     p2 = patch(vertex(ii,1), vertex(ii,2),phi(ii,1),'FaceColor','interp');
    fill3(vertex(ii,1), vertex(ii,2),phi(ii,1),phi(ii,1)); hold on;
end



