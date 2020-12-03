clear;

element = load('HW18.element');
vertex = load('HW18.vertex');

% for row=1:size(element,1)
%     ii = [element(row,1:3) element(row,1)];
%     p = patch(vertex(ii,1), vertex(ii,2), 'w', 'LineStyle','-');
% end

n1 = 416; %1
n2 = 417; %0

b = zeros(size(vertex,1),1);
A = sparse(size(vertex,1),size(vertex,1));



for ii = 1:size(element,1)
    
    
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
            A(element(ii,1),element(ii,1)) = A(element(ii,1),element(ii,1)) + (A2/L31) + ((-1)*A3/L12);
            A(element(ii,1),element(ii,2)) = A(element(ii,1),element(ii,2)) + (A3/L12);
            A(element(ii,1),element(ii,3)) = A(element(ii,1),element(ii,3)) + (-A2/L31);
            
            A(element(ii,2),element(ii,1)) = A(element(ii,2),element(ii,1)) + (A3/L12);
            A(element(ii,2),element(ii,2)) = A(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + ((-1)*A3/L12);
            A(element(ii,2),element(ii,3)) = A(element(ii,2),element(ii,3)) + (A1/L23);
            
            A(element(ii,3),element(ii,1)) = A(element(ii,3),element(ii,1)) + (-A2/L31);
            A(element(ii,3),element(ii,2)) = A(element(ii,3),element(ii,2)) + (A1/L23);
            A(element(ii,3),element(ii,3)) = A(element(ii,3),element(ii,3)) + (A2/L31) + ((-1)*A1/L23);
            
        elseif max(L) == L(1,1) %L12
            %             A3(-)
            A(element(ii,1),element(ii,1)) = A(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + (A3/L12);
            A(element(ii,1),element(ii,2)) = A(element(ii,1),element(ii,2)) + (-A3/L12);
            A(element(ii,1),element(ii,3)) = A(element(ii,1),element(ii,3)) + (A2/L31);
            
            A(element(ii,2),element(ii,1)) = A(element(ii,2),element(ii,1)) + (-A3/L12);
            A(element(ii,2),element(ii,2)) = A(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + (A3/L12);
            A(element(ii,2),element(ii,3)) = A(element(ii,2),element(ii,3)) + (A1/L23);
            
            A(element(ii,3),element(ii,1)) = A(element(ii,3),element(ii,1)) + (A2/L31);
            A(element(ii,3),element(ii,2)) = A(element(ii,3),element(ii,2)) + (A1/L23);
            A(element(ii,3),element(ii,3)) = A(element(ii,3),element(ii,3)) + ((-1)*A2/L31) + ((-1)*A1/L23);
            
        elseif max(L) == L(1,2) %L23
            %             A1(-)
            
            A(element(ii,1),element(ii,1)) = A(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + ((-1)*A3/L12);
            A(element(ii,1),element(ii,2)) = A(element(ii,1),element(ii,2)) + (A3/L12);
            A(element(ii,1),element(ii,3)) = A(element(ii,1),element(ii,3)) + (A2/L31);
            
            A(element(ii,2),element(ii,1)) = A(element(ii,2),element(ii,1)) + (A3/L12);
            A(element(ii,2),element(ii,2)) = A(element(ii,2),element(ii,2)) + (A1/L23) + (-A3/L12);
            A(element(ii,2),element(ii,3)) = A(element(ii,2),element(ii,3)) + (-A1/L23);
            
            A(element(ii,3),element(ii,1)) = A(element(ii,3),element(ii,1)) + (A2/L31);
            A(element(ii,3),element(ii,2)) = A(element(ii,3),element(ii,2)) + (-A1/L23);
            A(element(ii,3),element(ii,3)) = A(element(ii,3),element(ii,3)) + (-A2/L31) + (A1/L23);
            
        end
        
    else
        %���� �Ǵ� ���� �ﰢ���� ��
        
        A(element(ii,1),element(ii,1)) = A(element(ii,1),element(ii,1)) + ((-1)*A2/L31) + ((-1)*A3/L12);
        A(element(ii,1),element(ii,2)) = A(element(ii,1),element(ii,2)) + (A3/L12);
        A(element(ii,1),element(ii,3)) = A(element(ii,1),element(ii,3)) + (A2/L31);
        
        A(element(ii,2),element(ii,1)) = A(element(ii,2),element(ii,1)) + (A3/L12);
        A(element(ii,2),element(ii,2)) = A(element(ii,2),element(ii,2)) + ((-1)*A1/L23) + ((-1)*A3/L12);
        A(element(ii,2),element(ii,3)) = A(element(ii,2),element(ii,3)) + (A1/L23);
        
        A(element(ii,3),element(ii,1)) = A(element(ii,3),element(ii,1)) + (A2/L31);
        A(element(ii,3),element(ii,2)) = A(element(ii,3),element(ii,2)) + (A1/L23);
        A(element(ii,3),element(ii,3)) = A(element(ii,3),element(ii,3)) + ((-1)*A2/L31) + ((-1)*A1/L23);
    end
    
end


b(n1,1) = 1;
A(n1,:) =0;
A(n1,n1) =1;

b(n2,1) = 0;
A(n2,:) =0;
A(n2,n2) =1;

phi = A \ b;



for row=1:size(element,1)
    ii = [element(row,1:3)];      
    fill3(vertex(ii,1), vertex(ii,2),phi(ii,1),phi(ii,1)); hold on;
end

