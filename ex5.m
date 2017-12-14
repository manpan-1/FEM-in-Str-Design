% Exercise 4: Frame solved using corotational elements with plasticity.

% Set globals
clear all;
global wx wy E Et yield H xg yg  

%% weight factors
wx=[0.5 0.5];       xg=[0.21132 0.78868];
wy=[0.12948 0.27971 0.38183 0.41796 0.38183 0.27971 0.12948];    yg=[-0.94911 -0.74153 -0.40585 0 0.40585 0.74153 0.94911];

% Initialise variables for the output
forc = [];
dispx = [];
dispy = [];

% Define approximation tolerance for the residual
tol = 1E-4; 

% Cross-section
L = 120;
a = 3;
b = 2;

% Material properties
E = 720;
Et = E/10;
yield = 10.44;
H = Et/(1-Et/E);

% Structural properties
I = a*b^3/12;
A = a*b;
EA = E*A;
EI = E*I;

% Geometrical parameters of the structure (elements, nodes etc.)
coor = zeros(11,8);
elem = zeros(10,2);

Ne = size(elem,1);
Np = size(coor,1);

%       x   y   BC_x  BC_y    BC_xy  Load_x Load_y M_xy
coor = [0   0   1     1       0      0      0      0    
        0   24  0     0       0      0      0      0
        0   48  0     0       0      0      0      0
        0   72  0     0       0      0      0      0
        0   96  0     0       0      0      0      0
        0   120 0     0       0      0      0      0
        24  120 0     0       0      0     -1e7    0
        48  120 0     0       0      0      0      0
        72  120 0     0       0      0      0      0
        96  120 0     0       0      0      0      0
        120 120 1     1       0      0      0      0];
  
elem = [1, 2
        2, 3
        3, 4
        4, 5
        5, 6
        6, 7
        7, 8 
        8, 9
        9, 10
        10, 11];

afg = [];
P = [];

% Build an array of the active degrees of freedom (acc. to boundary 
% conditions)
for i = 1:Np
    if coor(i, 3) == 0 
       afg = [afg; 3*(i - 1) + 1];
       P = [P; coor(i, 6)];       
    end
    if coor(i, 4) == 0 
       afg = [afg; 3*(i - 1) + 2];
       P = [P; coor(i, 7)];       
    end    
    if coor(i, 5) == 0 
       afg = [afg; 3*(i - 1) + 3];
       P = [P; coor(i, 8)];       
    end    
end  

nfg = size(afg,1);

%% assembling 
% construct global matrices
Kg = zeros(3*Np, 3*Np); % size of matrix determined by number of degrees of freedom per node
qg = zeros(3*Np, 1); % internal force vector
d = zeros(3*Np, 1); %displacement matrix
do = d; %previous step displacement (initially zero)
Pr = d; %full matrix of external forces
residualg = d; %residual forces global matrix


sgo=zeros(7*Ne,2); %stress matrix for gausspoints
epsgo=zeros(7*Ne,2);
ego=zeros(7*Ne,2);
sn=zeros(7,2);
epsn=zeros(7,2);


% initialisation
for i = 1:Ne  
  m1 = elem(i, 1);
  m2 = elem(i, 2);
  
  % coordinates matrix
  x = [coor(m1,1:2), coor(m2,1:2)]';
  v = [3*m1-2:3*m1, 3*m2-2:3*m2]';
  gp = [7*m1-6:7*m2-7]';
  [q, K, sn, epsn] = corotbeamplastic ( a, b, x, d(v), sgo(gp,1:2) ,epsgo(gp,1:2), do(v) );
  Kg(v, v) = Kg(v, v) + K;
  
end

% Reduced (for BCs) global tangent stiffness matrix.
Kr = Kg(afg, afg);  %initial tangent stiffness matrix.

% displacement increment.
d_inc = -0.1;

% Auxiliary vector for the displacement control.
aux = [zeros(17, 1); 1 ; zeros(11, 1)];
dd = d_inc*aux;

%% apply displacement control 
% define tangent stiffnes matrix
Ka = [Kr, -P; aux', 0];

% lambda
la = 0;
max_idx = 18;

% perform analysis
step3 = 0;

while d(19) < 100
    % assemble tangent vector 
    t = [Kr\P; 1];
    
    % predictor
    %identify maximum displacement: node number and direction
    if d(19) < 90 && abs(d(20)) > 48
        max_idx=17;
        d_inc=0.1;
    end 
    
    % Change displacement control to the node with the highest displacement,
    aux = [zeros(max_idx-1, 1); 1 ; zeros(length(t) - max_idx - 1, 1) ]; 
    
    dd = aux*d_inc;

    predictor = dd(max_idx)/t(max_idx)*t;
    
    % Write dlambda and ddisplacement
    dd = predictor(1:end-1);
    dF = predictor(end);
      
    % Update displacement vector
    d(afg, 1) = d(afg, 1)+dd;
    
    % Update lambda
    la = la+dF;
    
    % Create the big matrix
    % Pr(afg,1) = Pr(afg,1)+P;
    r = 100;
    Ka = [Kr, -P; aux', 0];
    
    while r > tol
        %assemble
        Kg = zeros(3*Np, 3*Np);
        qg = zeros(3*Np, 1);
        sgn=zeros(7*Ne,2);
        epsgn=zeros(7*Ne,2);
        for i = 1:Ne
            m1 = elem(i, 1);
            m2 = elem(i, 2);
            
            %coordinates matrix
            x = [coor(m1, 1:2), coor(m2, 1:2)]';
            v = [3*m1 - 2:3*m1, 3*m2 - 2:3*m2]';     
            
            % Local stiffness, load vector and strains
            [q, K, sn, epsn] = corotbeamplastic( a, b, x, d(v),sgo(gp,1:2), epsgo(gp,1:2), do(v));
            
            % Assemble to global
            Kg(v, v) = Kg(v, v) + K;
            qg(v, 1) = qg(v, 1) + q;
            sgn(gp,1:2)=sgn(gp,1:2)+sn;
            epsgn(gp,1:2)=epsgn(gp,1:2)+epsn;
        end
        
        qr = qg(afg, 1);
        Kr = Kg(afg, afg);  %converged stiffness matrix
        
        residualr = qr - la*P;
        r = norm(residualr);
        if r > tol
            um  =  -Ka\[residualr; 0];
            dd = um(1:end - 1);   
            dF = um(end);
            d(afg, 1) = d(afg, 1) + dd;
            la = la + dF;
        end
    end
    
    % store matrices at the current converged step to use on the next step.
    do= d;
    sgo= sgn;
    epsgo= epsgn;
    
    % keep values of force and displacement at the control node for later
    % plotting.
    forc = [forc, la];
    dispx = [dispx, d(19)];
    dispy = [dispy, d(20)];
end
% 
% P = 1;
% u = d(19);
% v = -d(20);
% 
% %defbeam(coor, elem, d, 0.1)
% 
% figure;
% plot(dispx, forc, -dispy, forc);

defbeam(coor,elem,d,1)
