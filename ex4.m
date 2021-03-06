% Exercise 4: Frame solved using corotational elements.

% Initialise variables for the output
forc = [];
dispx = [];
dispy = [];

% Define approximation tolerance for the residual
tol = 1E-3; 

% Cross-section
L = 120;
a = 2;
b = 3;

% Structural properties
E = 720;
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
        24  120 0     0       0      0     -3      0
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
Pr = d; %full matrix of external forces
residualg = d; %residual forces global matrix

% initialisation
for i = 1:Ne  
  m1 = elem(i, 1);
  m2 = elem(i, 2);
  
  % coordinates matrix
  x = [coor(m1,1:2), coor(m2,1:2)]';
  v = [3*m1-2:3*m1, 3*m2-2:3*m2]';  
  [q, K] = corotbeam(EA, EI, x, d(v));
  Kg(v, v) = Kg(v, v) + K;
  
end

% Reduced (for BCs) global tangent stiffness matrix.
Kr = Kg(afg, afg);  %initial tangent stiffness matrix.

% displacement increment.
d_inc = -0.05;

% Auxiliary vector for the displacement control.
aux = [zeros(17, 1); 1 ; zeros(11, 1)];
dd = d_inc*aux;

%% apply displacement control 
% define tangent stiffnes matrix
Ka = [Kr, -P; aux', 0];

% lambda
la = 0;
max_idx = 18;

% direction vector
s = 1;

% perform analysis
step3 = 0;

while la < 2.0
    % assemble tangent vector 
    t = [Kr\P; 1];
    
    % predictor
    if abs(d(20)) < 60 && max_idx ~= 17  
        max_idx = 18;               
    elseif d(19) < 93 && step3 == 0
        max_idx = 17;
        d_inc = 0.02;
    end 
    
    if abs(d(19)) > 93    
        max_idx = 21;
        d_inc = -0.1;
        step3 = 1;
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
        for i = 1:Ne
            m1 = elem(i, 1);
            m2 = elem(i, 2);
            
            %coordinates matrix
            x = [coor(m1, 1:2), coor(m2, 1:2)]';
            v = [3*m1 - 2:3*m1, 3*m2 - 2:3*m2]';                   
            [q, K] = corotbeam(EA, EI, x, d(v));
            Kg(v, v) = Kg(v, v) + K;
            qg(v, 1) = qg(v, 1) + q;
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
    
    forc = [forc, la];
    dispx = [dispx, d(19)];
    dispy = [dispy, d(20)];
end

P = 1;
u = d(19);
v = -d(20);

defbeam(coor, elem, d, 0.1)

figure;
plot(dispx, forc, -dispy, forc);

