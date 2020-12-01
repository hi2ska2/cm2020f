% Draw.m

clear all;
close all;

ele = load('elements.txt');
ver = load('vertex.txt');
length = zeros(3,1);
xdif = zeros(3,1);
ydif = zeros(3,1);
cnall = zeros(size(ele,1),3);
figure(1);
phi = zeros(size(ver,1),1);
B = zeros(size(ver,1), size(ver,1));
L = zeros(size(ver,1), size(ver,1));
C = zeros(size(ver,1), size(ver,1));
linked_elements = zeros(size(ele,1),14);
count2 = 0;
b = zeros(size(ver,1),1);

% 한 Vertex 기준으로 연결 된 vertices 의 최대 값
for row = 1:size(ele, 1)
    
   
    count = 0;
    for row2 = 1:size(ele,1)
        x1 = ele(row2, 1);
        x2 = ele(row2, 2);
        x3 = ele(row2, 3);
        if x1 == row
            
            count = count+1;
            linked_elements(row, count) = row2;
        elseif x2== row
            count = count+1;
            linked_elements(row, count) = row2;
        elseif x3 == row
            count = count+1;
            linked_elements(row, count) = row2;
        end
    end
    if count > count2
    count2 = count;
    else
    end
end
sort_ele = zeros(1564,3);




  
l_e = linked_elements;
number = zeros(size(ele,1),1); % nth vertex 를 둘러 싸고 있는 삼각형의 갯수
for row = 1:size(ver,1)
    for i = 1:14
        if l_e(row, i) == 0
            %l_e(row,i) = linked_elements(row,1);
            number(row,1) = i-1;
            break;
        end
        
    end
end






for row = 1:size(ele, 1)
    n = number(row,1);
    ii = [ele(row, 1:3) ele(row,1)];
    jj = [l_e(row, 1:n)];
    %hold on;
    %p1 = patch(ver(ii,1), ver(ii,2), 'r', 'LineStyle', '-');
    %p2 = patch(cnall(jj,1), cnall(jj,2), 'w', 'LineStyle', '-');
    %hold on;
    %{
    xdif(1,1) = ver(ii(1,1),1)-ver(ii(1,2),1);
    xdif(2,1) = ver(ii(1,2),1)-ver(ii(1,3),1);
    xdif(3,1) = ver(ii(1,3),1)-ver(ii(1,1),1);
    ydif(1,1) = ver(ii(1,1),2)-ver(ii(1,2),2);
    ydif(2,1) = ver(ii(1,2),2)-ver(ii(1,3),2);
    ydif(3,1) = ver(ii(1,3),2)-ver(ii(1,1),2);
    length = sqrt(xdif.^2 + ydif.^2);
    cos = (length(1,1)^2 +length(2,1)^2 - length(3,1)^2)/(2*length(1,1)*length(2,1));
    theta = acos(cos);
    %}
    
    cor = [ ver(ii(1,1),1) ver(ii(1,2),1) ver(ii(1,3),1)  ; ver(ii(1,1),2) ver(ii(1,2),2) ver(ii(1,3),2) ];
    [r, cn]=circumcircle(cor, 0);
    %plot(cn(1,1),cn(2,1), '.');
    cnall(row,1) = cn(1,1);
    cnall(row,2) = cn(2,1);
    
end


%%%%% l_e sorting %%%%%
l_e_update = zeros(size(ele,1),14);
baby=zeros(2,1);
angle = zeros(size(ver,1),13);
for row = 1:size(ver,1)
    centroid = zeros(13,2);
    vector = zeros(13,2);
    n = number(row,1);
    %angle = zeros(13,1);
    angle_sorted = zeros(2,1);
    stvector = zeros(1,2);
    
    for i = 1:n
        
        centroid(i,1) = (ver(ele(l_e(row,i),1),1) + ver(ele(l_e(row,i),2),1) + ver(ele(l_e(row,i),3),1))/3; %centroid x position 
        centroid(i,2) = (ver(ele(l_e(row,i),1),2) + ver(ele(l_e(row,i),2),2) + ver(ele(l_e(row,i),3),2))/3; %centroid y position
        vector(i,1) = centroid(i,1) - ver(row,1);
        vector(i,2) = centroid(i,2) - ver(row,2);
        stvector(1,1) = 1e-9;
        costheta = stvector(1,1)*vector(i,1)/(norm(stvector)*norm(vector(i,:)));
        if vector(i,2) > 0
            angle(row,i) = acos(costheta);
        else
            angle(row, i) = 2*pi - acos(costheta);
        end
    end
    
    for i = 1:n
        for k = 1:(n-i)
            if angle(row,i) > angle(row,i+k)
                angle_sorted(1,1) = angle(row,i);
                angle_sorted(2,1) = angle(row,i+k);
                angle(row,i) = angle_sorted(2,1);
                angle(row, i+k) = angle_sorted(1,1);
                baby(1,1) = l_e(row,i);
                baby(2,1) = l_e(row, i+k);
                l_e(row,i) = baby(2,1);
                l_e(row,i+k) = baby(1,1);%element 재배정
            
            end
        end
    end
    

end
%sorting 2 (끊어진 부분을 양 끝으로 2점이 겹치지 않는 삼각형이 양 끝에 있어야 한다)

for row = 1:size(ver,1)
    
    n = number(row,1);
    for i =1:n
        if size(intersect(ele(l_e(row, n),:),ele(l_e(row, 1),:)),2) == 2 %이 경우 1번과 2번은 서로 맞닿아 있다는 이야기
            temp = l_e(row,1);
            for j = 1:(n-1)
            
            l_e(row,j)=l_e(row,j+1);
            end
            l_e(row,n) = temp;

        elseif size(intersect(ele(l_e(row, n),:),ele(l_e(row, 1),:)),2) == 1 %양 끝 삼각형이 서로 선분을 공유 하지 않으면 제대로 sorting 되어있다는 뜻
            break;
        end
    end
end








u1 = zeros(1,3);
v1 = zeros(1,3);
u2 = zeros(1,3);
v2 = zeros(1,3);

for row = 1:size(ver,1)  %%% L 과 A 행렬 구하기
    n = number(row,1);
    for j = 1:n
        if (ele(l_e(row,j),1) == row)
            i1 = ele(l_e(row,j),2);
            i2 = ele(l_e(row,j),3);
        elseif (ele(l_e(row,j),2) == row)
            i1 = ele(l_e(row,j),1);
            i2 = ele(l_e(row,j),3);
        elseif (ele(l_e(row,j),3) == row)
            i1 = ele(l_e(row,j),1);
            i2 = ele(l_e(row,j),2);
        end
        L(row, i1) = norm(ver(i1,:)-ver(row,:));
        L(row, i2) = norm(ver(i2,:)-ver(row,:));
        u1(1,:) = cnall(i1, :) - ver(row,:);
        v1(1,:) = ver(i1, :) - ver(row,:);
        u2(1,:) = cnall(i2, :) - ver(row,:);
        v2(1,:) = ver(i2, :) - ver(row,:);
        theta1 = acos((u1*transpose(v1))/(norm(u1)*norm(v1)));
        theta2 = acos((u2*transpose(v2))/(norm(u2)*norm(v2)));
        d1=abs(norm(u1)*sin(theta1));
        d2=abs(norm(u2)*sin(theta2));
        B(row, i1) = d1;%A_row,j 구한것
        B(row, i2) = d2;
    end
end

for row = 1:size(ver,1) 
    n = number(row,1);
    for j = 1:n
        if (ele(l_e(row,j),1) == row)
            i1 = ele(l_e(row,j),2);
            i2 = ele(l_e(row,j),3);
        elseif (ele(l_e(row,j),2) == row)
            i1 = ele(l_e(row,j),1);
            i2 = ele(l_e(row,j),3);
        elseif (ele(l_e(row,j),3) == row)
            i1 = ele(l_e(row,j),1);
            i2 = ele(l_e(row,j),2);
        end
        C(row, row) = C(row, row) - B(row, i1)/L(row,i1) - B(row, i2)/L(row,i2);
        C(row, i1) = C(row, i1) + B(row, i1)/L(row,i1);
        C(row, i2) = C(row, i2) + B(row, i2)/L(row,i2);
        
        
    end
end

%Boundary condition 주기%
C(1,:) = 0;
C(1,1) = 1;
C(100,:) = 0;
C(100,100) = 1;
b(1,1) = 1;

phi = (-C)\b;



for row = 1:size(ele, 1)
    
    plot3(ver(row,1), ver(row,2), phi(row,1), 'o');
    hold on;
end


%{
for row = 3:3%size(ele, 1)
    
    ii = [ele(row, :)];
    trisurf(ii, ver(ii,1), ver(ii,2),phi(ii,1));
    hold on;
    %triplot(ii, ver(ii,1), ver(ii,2));
end



%}
