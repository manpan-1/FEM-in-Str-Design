% example 2
% Lee's frame 10 elements
% linear analysis

clear all;

global EA EI

tol=1E-5; 

% structure

L=120;
a=3;
b=2;

E=720;
I=a*b^3/12;
A=a*b;
EA=E*A;
EI=E*I;

coor=zeros(11,8);
elem=zeros(10,2);

Ne=size(elem,1);
Np=size(coor,1);

%     x   y  BC_x  BC_y BC_xy Load_x Load_y Load_xy
coor=[0   0   1 1 0 0  0   0    
      0   24  0 0 0 0  0   0
      0   48  0 0 0 0  0   0
      0   72  0 0 0 0  0   0
      0    96 0 0 0 0  0   0
      0   120 0 0 0 0  0   0
      24  120 0 0 0 0 -1   0
      48  120 0 0 0 0  0   0
      72  120 0 0 0 0  0   0
      96  120 0 0 0 0  0   0
      120 120 1 1 0 0  0   0];
  
  elem=[1 2
        2 3
        3 4
        4 5
        5 6
        6 7
        7 8 
        8 9
        9 10
        10 11];

afg=[];
P=[];

for i=1:Np
    if coor(i,3)==0 
       afg=[afg;3*(i-1)+1];
       P=[P;coor(i,6)];       
    end
    if coor(i,4)==0 
       afg=[afg;3*(i-1)+2];
       P=[P;coor(i,7)];       
    end    
    if coor(i,5)==0 
       afg=[afg;3*(i-1)+3];
       P=[P;coor(i,8)];       
    end    
end  

nfg=size(afg,1);

% assembling

Kg=zeros(3*Np,3*Np);

for i=1:Ne
  
  m1=elem(i,1);
  m2=elem(i,2);
  x=[coor(m1,1:2) coor(m2,1:2)]';
  v=[3*m1-2:3*m1 3*m2-2:3*m2]';
 
  L=sqrt((x(3)-x(1))^2+(x(4)-x(2))^2);
  
  c=(x(3)-x(1))/L;
  s=(x(4)-x(2))/L;
  %local stiffness matrix
  ke=[EA/L   0           0         -EA/L   0           0
    	 0      12*EI/L^3   6*EI/L^2   0     -12*EI/L^3   6*EI/L^2
    	 0      6*EI/L^2    4*EI/L     0     -6*EI/L^2    2*EI/L
   	 -EA/L   0           0          EA/L   0           0
    	 0     -12*EI/L^3  -6*EI/L^2   0      12*EI/L^3  -6*EI/L^2
    	 0      6*EI/L^2    2*EI/L     0     -6*EI/L^2    4*EI/L];
  %transformation matrix
  T=[c  s  0  0  0  0
  	-s  c  0  0  0  0
   	 0  0  1  0  0  0
   	 0  0  0  c  s  0
   	 0  0  0 -s  c  0
   	 0  0  0  0  0  1];
  %global coordinates stiffness matrix 
  ke1=T'*ke*T;
  %global stifness matrix
  Kg(v,v)=Kg(v,v)+ke1;  
  
end 

Kr=Kg(afg,afg);

dr=Kr\P;

d=zeros(3*Np);
d(afg)=dr;

P=1
u=d(19)
v=-d(20)


defbeam(coor,elem,d,1)

    
 