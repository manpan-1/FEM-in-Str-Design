% Exercise 6: Non-linear plastic analysis of thin ring using Green's strain.

%% Define globals
% The following globals are used by the pstress2d function
global E nu H yield

%% Input
% Define the approximation tolerance.
tol = 1E-4;

% Load
load = -100000;

% Define step, percentage for the given load.
% 0<step<1
d_step = -0.5;

%% Structure
% Define material properties.
E = 1E7;
yield = 2e5;
E_t = E/10;
H = E_t/(1-E_t/E);
nu = 0.25;
th = 1;

% Radius and XXX of the ring.
R = 20;
d = 1/2;

% Polar angle of nodes (csys origin at the circle centre).
an = pi/(2*40);

% Initialise matrices for xy-coordinates and element connectivity.
coor = zeros(41*3,6);
elem = zeros(40*2,4);

% Get the number of elements and the number of nodes.
Ne = size(elem,1);
Np = size(coor,1);

% Populate the nodal coordinate matrix.
for i = 1:41
    for j = 0:2
        ne = i+j*(41);
        ang = an*(i-1);
        coor(ne,1:2) = [(R+d*j)*cos(ang),(R+d*j)*sin(ang)];
    end
end

% Populate the element connectivity table.
for i = 1:40
    for j = 0:1
        ne = i+j*40;
        n1 = i+j*(41);
        elem(ne,:) = [n1 n1+41 n1+42 n1+1];
    end
end

% BC: fix x and y at the lower-right face of the ring.
coor(1,3:4) = [1 1];
coor(42,3:4) = [1 1];
coor(83,3:4) = [1 1];

% BC: fix x at the upper-left face of the ring.
coor(41,3) = 1;
coor(82,3) = 1;
coor(123,3) = 1;

% Apply vertical load on the uppermost node.
coor(123,6) = load;

% Plot the undeformed structure
figure;
hold on;
for i = 1:Ne
    n1 = elem(i,1);
    n2 = elem(i,2);
    n3 = elem(i,3);
    n4 = elem(i,4);
    x_i = [coor(n1,1),coor(n2,1),coor(n3,1),coor(n4,1),coor(n1,1)];
    y_i = [coor(n1,2),coor(n2,2),coor(n3,2),coor(n4,2),coor(n1,2)];
    plot(x_i, y_i, 'b');
end
grid;
axis equal;

% Initialise arrays to store the active degrees of freadom table and the
% load vector.
afg = [];
vload = [];

% Populate the active DOF's table and force vector according to the info
% given in the 'coor' table.
for i = 1:Np
    if coor(i,3) == 0
        afg = [afg;2*(i-1)+1];
        vload = [vload;coor(i,5)];
    end
    if coor(i,4) == 0
        afg = [afg;2*(i-1)+2];
        vload = [vload;coor(i,6)];
    end
end

%% Path following.
%
% Initialise the internal force vector, the tangent stiffness matrix and
% th displacement vector with zeros.
qt = zeros(2*Np, 1);
kt = zeros(2*Np, 2*Np);
d = zeros(2*Np, 1);

% Initialise stress and strain state tables and previous step displacements
%(all zero for the first step).
sgo = cell(Ne, 1);
sgn = cell(Ne, 1);
epsgo = zeros(Ne, 4);
epsgn = zeros(Ne, 4);
do = d;

% First prediction for internal force vector and tangent stiffness for zero
% displacements.
% Loop through the elements and assemble q and K matrices.
for i = 1:Ne
    % The 4 nodes of the current element.
    m1 = elem(i, 1);
    m2 = elem(i, 2);
    m3 = elem(i, 3);
    m4 = elem(i, 4);
    
    % xy-coordinates of the 4 nodes of the current element.
    xy = [coor(m1,1:2), coor(m2,1:2), coor(m3,1:2), coor(m4,1:2)];
    vu = [2*m1-1:2*m1,  2*m2-1:2*m2,  2*m3-1:2*m3,  2*m4-1:2*m4];
    
    % Element internal force and tangent stiffness matrix on the global csys.
    [q_e, k_e] = isoplnonlin(E, nu, th, xy, zeros(8, 1));
    
    % Add them to the global matrices.
    qt(vu) = qt(vu) + q_e;
    kt(vu,vu) = kt(vu,vu) + k_e;
    
    % Fill the sgo and epsgo matrices with zeros (for the initial step)
    sgo{i} = {zeros(3, 1), zeros(3, 1), zeros(3, 1), zeros(3, 1)};
    sgn{i} = {zeros(3, 1), zeros(3, 1), zeros(3, 1), zeros(3, 1)};
end

% Reduced tangent matrices(for the active DOF's).
qr = qt(afg);
kr = kt(afg,afg);

% Auxiliary vector.
h_v = [zeros(1, 236), 1];

% Extended tangent stiffness matrix for the displacement control.
K_a = [kr, -vload; h_v, 0];

% Initialise an array to store monitor the vertical displacement of the
% loaded node and the load factor, lambda.
disp_hist = zeros(ceil(1/d_step), 1);
lmbd_hist = zeros(ceil(1/d_step), 1);
la = 0;
i = 0;

% Path following loop.
while la<1;
    i = i+1;

    t = [kr\vload; 1];
    
    % Initial predictor and displacement vector.
    predictor = d_step/t(237)*t;
    dd = predictor(1:237);
    
    % Initial load prediction.
    dF = predictor(238);
    la = la + dF;
    lmbd_hist(i) = la;
        
    % Update displacements on the reduced matrix.
    d(afg) = d(afg)+dd;
    
    % Approximate the target load of the current step.
    r = 100;
    while r>tol;

        % Reset matrices
        kt = zeros(2*Np, 2*Np);
        qt = zeros(2*Np, 1);
        
        % Loop through the elements and assemble q and K matrices.
        for j = 1:Ne
            % The 4 nodes of the current element.
            m1 = elem(j, 1);
            m2 = elem(j, 2);
            m3 = elem(j, 3);
            m4 = elem(j, 4);
            
            % xy-coordinates and positions in the global matrix of the 4
            % nodes of the current element.
            xy = [coor(m1,1:2), coor(m2,1:2), coor(m3,1:2), coor(m4,1:2)];
            vu = [2*m1-1:2*m1,  2*m2-1:2*m2,  2*m3-1:2*m3,  2*m4-1:2*m4];
            
            % Element internal force and tangent stiffness matrix on the global csys.
            % The current step sgo and epsgo get updated.
            [q_e, k_e, sgn{j}, epsgn(j, :)] = isoplnonlinplast(th, xy, d(vu), sgo{j}, epsgo(j, :), do(vu));
            
            % Add them to the global matrices.
            qt(vu) = qt(vu) + q_e;
            kt(vu,vu) = kt(vu,vu) + k_e;
        end
        
        % Reduced matrices.
        qr = qt(afg);
        kr = kt(afg,afg);
        
        % Update extended tangent stiffness matrix for displacement
        % control.
        K_a = [kr, -vload; h_v, 0];
        
        % Update residual.
        resid = qr - la*vload;
        r = norm(resid);
                
        % Update the predictor and displacements.
        predictor = -K_a\[resid;0];
        dd = predictor(1:237);
        d(afg) = d(afg)+dd;
        
        % Update load step.
        dF = predictor(238);
        la = la + dF;
        
    end
    disp_hist(i) = d(82);
    
    % Update stored converged data (d, sg, eps)for plasticity.
    do = d;
    sgo = sgn;
    epsgo = epsgn;
end

% Global displacements.
d(afg) = d(afg)+dd;

%% Post-rocessing
% Print on screen the requested displacement.
disp(d(82))

% Define scale factor for printing the deformed shape.
sc = 0.1;

% Plot the deformed shape.
for i = 1:Ne
    n1 = elem(i,1);
    n2 = elem(i,2);
    n3 = elem(i,3);
    n4 = elem(i,4);
    xy = [coor(n1,1)+d(2*n1-1)*sc,coor(n2,1)+d(2*n2-1)*sc,coor(n3,1)+d(2*n3-1)*sc,coor(n4,1)+d(2*n4-1)*sc,coor(n1,1)+d(2*n1-1)*sc];
    y = [coor(n1,2)+d(2*n1-0)*sc,coor(n2,2)+d(2*n2-0)*sc,coor(n3,2)+d(2*n3-0)*sc,coor(n4,2)+d(2*n4-0)*sc,coor(n1,2)+d(2*n1-0)*sc];
    plot(xy,y,'r');
end

% Release figure
hold off;

% Plot load displacement curve
figure;
plot([0, -disp_hist], -[0, lmbd_hist]*load, '-o')

