% example 3
% thin ring - linear analysis

clear all;

E=1E7;
nu=0.25;
th=1;

R=20;
d=1/2;

an=pi/2/40;

coor=zeros(41*3,6);
elem=zeros(40*2,4);

Ne=size(elem,1);
Np=size(coor,1);

for i=1:41
  for j=0:2  
    ne=i+j*(41);  
    ang=an*(i-1);
    coor(ne,1:2)=[(R+d*j)*cos(ang),(R+d*j)*sin(ang)];
  end
end  

coor(1,3:4)=[1 1];
coor(42,3:4)=[1 1];
coor(83,3:4)=[1 1];

coor(41,3)=1;
coor(82,3)=1;
coor(123,3)=1;

coor(123,6)=-5000;

for i=1:40
  for j=0:1
     ne=i+j*40;
     n1=i+j*(41);
     elem(ne,:)=[n1 n1+41 n1+42 n1+1];
  end
end  

hold on;

for i=1:Ne
  n1=elem(i,1);
  n2=elem(i,2);
  n3=elem(i,3);
  n4=elem(i,4);
  x=[coor(n1,1),coor(n2,1),coor(n3,1),coor(n4,1),coor(n1,1)];
  y=[coor(n1,2),coor(n2,2),coor(n3,2),coor(n4,2),coor(n1,2)];
  plot(x,y,'b');
end
grid;
axis equal;

afg=[];
vload=[];

for i=1:Np
    if coor(i,3)==0 
       afg=[afg;2*(i-1)+1];
       vload=[vload;coor(i,5)];       
    end
    if coor(i,4)==0 
       afg=[afg;2*(i-1)+2];
       vload=[vload;coor(i,6)];       
    end    
end  

kt=zeros(2*Np,2*Np);

for i=1:Ne
  
  m1=elem(i,1);
  m2=elem(i,2);
  m3=elem(i,3);
  m4=elem(i,4);
 
  x=[coor(m1,1:2) coor(m2,1:2) coor(m3,1:2) coor(m4,1:2)];
  v=[2*m1-1:2*m1 2*m2-1:2*m2 2*m3-1:2*m3 2*m4-1:2*m4];
  
  [k]=isopl(E,nu,th,x);
  
  kt(v,v)=kt(v,v)+k;  
  
end    

kr=kt(afg,afg);

dr=kr\vload;

d=zeros(2*Np,1);
d(afg)=dr;

d(82)

sc=1;

for i=1:Ne
  n1=elem(i,1);
  n2=elem(i,2);
  n3=elem(i,3);
  n4=elem(i,4);
  x=[coor(n1,1)+d(2*n1-1)*sc,coor(n2,1)+d(2*n2-1)*sc,coor(n3,1)+d(2*n3-1)*sc,coor(n4,1)+d(2*n4-1)*sc,coor(n1,1)+d(2*n1-1)*sc];
  y=[coor(n1,2)+d(2*n1-0)*sc,coor(n2,2)+d(2*n2-0)*sc,coor(n3,2)+d(2*n3-0)*sc,coor(n4,2)+d(2*n4-0)*sc,coor(n1,2)+d(2*n1-0)*sc];
  plot(x,y,'r');
end


