% this function plot the eigenmodes 
% for homework 9 a

% coor : coordinates x and y of the nodes
% elem : node 1 and node 2 of the elements
% d  : displacements and rotations at the nodes

function [rien]=plotmode(coor,elem,d);

Ne=size(elem,1);

hold on;

for i=1:Ne

   X=[]; 
   Y=[];
   
   n1=elem(i,1);
   n2=elem(i,2);
   
   x1=coor(n1,1);
   y1=coor(n1,2);
   x2=coor(n2,1);
   y2=coor(n2,2);
   
   lo=sqrt((x2-x1)^2+(y2-y1)^2);
   
   c=(x2-x1)/lo;
   s=(y2-y1)/lo;
   
   T=[c s 0  0 0 0 
     -s c 0  0 0 0
      0 0 1  0 0 0
      0 0 0  c s 0
      0 0 0 -s c 0
      0 0 0  0 0 1];
   
   TT=[c s
      -s c];
  
   v=[3*n1-2:3*n1 3*n2-2:3*n2];

   p=T*d(v);
   
   for j=0:100
      x=j/100*lo;
      f1=1-3*(x/lo)^2+2*(x/lo)^3;
      f2=x*(1-x/lo)^2;
      f3=3*(x/lo)^2-2*(x/lo)^3;
      f4=x^2/lo*(x/lo-1);
      u=(1-x/lo)*p(1)+x/lo*p(4);
      v=f1*p(2)+f2*p(3)+f3*p(5)+f4*p(6);
      U=TT'*[x+u;v];
      ug=U(1);
      vg=U(2);
      X=[X,x1+ug];
      Y=[Y,y1+vg];
   end
   
   plot(X,Y,'b');
   
end   

rien=0;

