% plot the deformation figure of a 2D beam structure
% in non linear analysis assuming a corotational
% description.

% coor : coordinates x and y of the nodes
% elem : node 1 and node 2 of the elements
% d  : displacements and rotations at the nodes

function [rien]=defbeam(coor,elem,d,scale);

d=d*scale;

Np=size(coor,1);
Ne=size(elem,1);


xv=[];
yv=[];
for i=1:Np
   xv=[xv,coor(i,1)];
   yv=[yv,coor(i,2)];
end      
plot(xv,yv,'r');

hold on;

for i=1:Ne

   n1=elem(i,1);
   n2=elem(i,2);
   
   x1=coor(n1,1);
   y1=coor(n1,2);
   x2=coor(n2,1);
   y2=coor(n2,2);
   
   v=[3*n1-2:3*n1 3*n2-2:3*n2];

   p=d(v);
         
   x21=[x2-x1
        y2-y1];  
     
   d21=[p(4)-p(1)
        p(5)-p(2)];
   
   lo=sqrt(x21'*x21);
   ln=sqrt((x21+d21)'*(x21+d21));
   
   cosbe=(x21(1)+d21(1))/ln;
   sinbe=(x21(2)+d21(2))/ln;
   
   sinal=(x21(1)*d21(2)-x21(2)*d21(1))/ln/lo;
   cosal=(x21'*(x21+d21))/ln/lo;

   if sinal>=0 & cosal>=0
      al=asin(sinal);
   end
   if sinal>=0 & cosal<0
      al=acos(cosal);
   end
   if sinal<0 & cosal>=0
      al=asin(sinal);
   end
   if sinal<0 & cosal<0
      al=-acos(cosal);
   end
   
   t1=p(3)-al;
   t2=p(6)-al;
   
   xv=[];
   yv=[];
   for j=0:100
      x=j/100*lo;
      xl=j/100*ln;
      yl=x*(1-x/lo)^2*t1+x^2/lo*(x/lo-1)*t2;
      xg=x1+p(1)+cosbe*xl-sinbe*yl;
      yg=y1+p(2)+sinbe*xl+cosbe*yl;
      xv=[xv,xg];
      yv=[yv,yg];
   end
   
    plot(xv,yv,'b');
   
end   

rien=0;

hold off
