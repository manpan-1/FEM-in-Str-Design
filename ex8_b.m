
clear all;

global EA EI

forc=[];
dispx=[];
dispy=[];

tol=1E-4; 

% structure

L=2;
a=0.1;
b=0.1;

E=2*10^11;
I=a*b^3/12;
A=a*b;
EA=E*A;
EI=E*I;

Pcr=pi^2*EI/L^2;
coor=zeros(11,8);
elem=zeros(10,2);

Ne=size(elem,1);
Np=size(coor,1);

%     x   y  BC_x  BC_y   BC_xy  Load_x Load_y M_xy
coor=[0      0   1     1       0      0      0     0    
      0.2    0   0     0       0      0      0     0
      0.4    0   0     0       0      0      0     0
      0.6    0   0     0       0      0      0     0
      0.8    0   0     0       0      0      0     0
      1.0    0   0     0       0      0   Pcr/500  0 
      1.2    0   0     0       0      0      0     0
      1.4    0   0     0       0      0      0     0
      1.6    0   0     0       0      0      0     0
      1.8    0   0     0       0      0      0     0
      2.0    0   0     1       0    -Pcr     0     0];
  
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

%% assembling 
%construct stiffness matrix
Kg=zeros(3*Np,3*Np); % size of matrix determined by number of degrees of freedom per node
qg=zeros(3*Np,1); % internal force vector
d=zeros(3*Np,1); %displacement matrix
%initialisation
for i=1:Ne  
  m1=elem(i,1);
  m2=elem(i,2);
  %coordinates matrix
  x=[coor(m1,1:2) coor(m2,1:2)]';
  v=[3*m1-2:3*m1 3*m2-2:3*m2]';  
  [q,K]=corotbeam ( EA, EI, x, d(v) );
  Kg(v,v)=Kg(v,v)+K;
  qr=qg(afg,1);
end
Kr=Kg(afg,afg);  %initial tangent stiffness matrix

%% step 1
%apply imperfection

%force control 
%lambda
la=0;
%load step 
dF = 0.01;
%stiffness matrix

a=[];
while la < 1.3
   
    residualr=qr-la*P;
    
    r=norm(residualr);
    
    while r>tol
        dd=-Kr\residualr;
        d(afg,1)=d(afg,1)+dd;
        
        %assemble        
        Kg=zeros(3*Np,3*Np);
        qg=zeros(3*Np,1);
        for i=1:Ne
            m1=elem(i,1);
            m2=elem(i,2);
            %coordinates matrix
            x=[coor(m1,1:2) coor(m2,1:2)]';
            v=[3*m1-2:3*m1 3*m2-2:3*m2]';                   
            [q,K]=corotbeam ( EA, EI, x, d(v) );
            Kg(v,v)=Kg(v,v)+K;
            qg(v,1)=qg(v,1)+q;
        end
        qr=qg(afg,1);
        Kr=Kg(afg,afg);  %converged stiffness matrix
        
        residualr=qr-la*P;
        r=norm(residualr);
        
    end
    forc=[forc,la];
    dispx=[dispx,d(31)];
    dispy=[dispy,d(17)];
    la=la+dF;
   
end

figure;
plot(dispy,forc)
figure;
scale = 1
defbeam(coor,elem,scale*d,1.0)