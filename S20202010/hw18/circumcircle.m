function [r,cn]=circumcircle(cor,pl)
%  circumcircle: Compute the circum radius, circum center of a  triangle and plot it
%  Definition: Circumcircle of a triangle is a circle going through the 
%  three of its vertices. It has a center at the point where 
%  the  perpendicular bisectors  of the triangle meets.
%  <Synopsis>
%    [r,cn] = circumcircle(cor)
%    [r,cn] = circumcircle(cor,pl)
%  <Input Parameters>
%    cor    --> coordinates of the three vertices of a triangle
%    given as cor=[x_cooridnates;y_cooridates]
%    x_cooridnates= 1 by 3 vector 
%    y_cooridnates= 1 by 3 vector 
%    pl --> plotting option, If pl=1, plots the triangle and circle, else 
%    it does not plot
%  <Output parameters>
%    r        --> radius of the circum circle
%    cn       --> center of the circumcircle 
%    Written by Bishnu P. Lamichhane, Aston University, UK
%    blamichha@yahoo.com
%    Last revised: 5th Nov, 2007
%-----------------------------------------------------------------------
  if (nargin==1) pl=0,end 
  %check the dimension of cor
  if (size(cor) ~= [2,3]) 
    error('Needs three vertices of a triangle');
  end
  %error  check computing the area
  ar=polyarea(cor(1,:),cor(2,:));
  if (ar<1e-30)  
    error('Degenerate triangle','Three points are almost collinear');
  end
  %compute the length of sides (AB, BC and CA) of the triangle
  c=norm(cor(:,1)-cor(:,2));
  a=norm(cor(:,2)-cor(:,3));
  b=norm(cor(:,1)-cor(:,3));
  %use formula: R=abc/(4*area) to compute the circum radius
  r=a*b*c/(4*ar);
  %compute the barycentric coordinates of the circum center
  bar=[a^2*(-a^2+b^2+c^2),b^2*(a^2-b^2+c^2),c^2*(a^2+b^2-c^2)];
  %convert to the real coordinates
  cn=bar(1)*cor(:,1)+bar(2)*cor(:,2)+bar(3)*cor(:,3);
  cn=cn/sum(bar);
  %plot the result
  if (pl==1) 
    th=linspace(0,2*pi);
    x=cn(1)+r*cos(th);
    y=cn(2)+r*sin(th);
    pn=[cor,cor(:,1)];plot(x,y,pn(1,:),pn(2,:),'-*');
  end

