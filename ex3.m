% Exercise 3: Example solved with rotated engineering strain and
% displacement control

% Define the approximation tolerance.
tol = 1E-4;

% Structure geometry.
d1 = 1;
d2 = 1;
h = 0.1;

% Nodal coordinates.
xy1 = [0 0 d1 h];
xy2 = [d1 h d1+d2 0];

% Initial cross-sectional area
Ao = 0.01;

% Young's modulus
E = 210000;

% Load vector.
F = [0; -1];

% initialisation
d = zeros(2,1);
d1 = zeros(4,1);
d2 = zeros(4,1);
la = 0;
forc = 0;
dispx = 0;
dispy = 0;

% Element formulation using rotational engineering strain
[fe1, ke1] = trussroteng(xy1, zeros(1, 4), E, Ao);
[fe2, ke2] = trussroteng(xy2, zeros(1, 4), E, Ao); 

%define internal loads matrix
feg = fe1(3:4)+fe2(1:2);
%define stiffnes matrix
keg=ke1(3:4,3:4)+ke2(1:2,1:2); 
%define tangent stiffnes matrix
Ka=[keg,-F;0,1,0];
dv2=-0.01;

for i=1:25
  %assemble tangent vector
  t=[keg\F;1];
  %predictor
  predictor=dv2/t(2)*t;
  dd=predictor(1:2);
  %update displacement
  d=d+dd;
  dF=predictor(3);
  la=la+dF;
  d1=[0 0 d(1) d(2)];
  d2=[d(1) d(2) 0 0];
  r=100;
  while r>tol
     %assemble the stiffness and loads matrix       
     [fe1,ke1]=trussroteng(xy1,d1,E,Ao);
     [fe2,ke2]=trussroteng(xy2,d2,E,Ao); 
     feg=fe1(3:4)+fe2(1:2);
     keg=ke1(3:4,3:4)+ke2(1:2,1:2); 
    %calculate residual
     resid=feg-la*F;
     r=norm(resid);
     if r>tol         
         %predictor based on initial tangent stiffness matrix
         predictor=-Ka\[resid;0];
         %update displacement         
         dd=predictor(1:2);
         d=d+dd;
         dF=predictor(3);
         la=la+dF;
         d1=[0 0 d(1) d(2)];
         d2=[d(1) d(2) 0 0];
     end
  end  
  
  forc=[forc,la];
  dispx=[dispx,d(1)];
  dispy=[dispy,d(2)];
  
end    


figure(1);
plot(dispx,forc,'r-x');
hold on;
grid;
figure(2);
plot(-dispy,forc,'r-x');
hold on;
grid;