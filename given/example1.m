% example1
% example page 3.2 with load control


tol=1E-4;

d1=1;
d2=0.5;
h=0.1;

Ao=0.01;
E=200000;

x1=[0 0 d1 h];
x2=[d1 h d1+d2 0];

F=[0
  -1];

% initialisation

d=zeros(2,1);
d1=zeros(4,1);
d2=zeros(4,1);
la=0;   

forc=[0];
dispx=[0];
dispy=[0];

[fe1,ke1]=trussgreen(x1,zeros(1,4),E,Ao);
[fe2,ke2]=trussgreen(x2,zeros(1,4),E,Ao); 
feg=fe1(3:4)+fe2(1:2);
keg=ke1(3:4,3:4)+ke2(1:2,1:2); 


dF=0.2;

for i=1:15
  
  i  
    
  la=la+dF;
  resid=feg-la*F;
  r=norm(resid);
 
  while r>tol;
      
     dd=-keg\resid;
     d=d+dd;
     
     d1=[0 0 d(1) d(2)];
     d2=[d(1) d(2) 0 0];
     
     [fe1,ke1]=trussgreen(x1,d1,E,Ao);
     [fe2,ke2]=trussgreen(x2,d2,E,Ao); 
     feg=fe1(3:4)+fe2(1:2);
     keg=ke1(3:4,3:4)+ke2(1:2,1:2); 
  
     resid=feg-la*F;
     r=norm(resid)
  end  
  
  forc=[forc,la];
  dispx=[dispx,d(1)];
  dispy=[dispy,d(2)];
  
end    

hold on;
figure(1);
plot(dispx,forc,'b-o');
grid;
figure(2);
plot(-dispy,forc,'b-o');
grid;