% Exercise 6: Non-linear analysis of thin ring using Green's strain.

%% Input
% Define the approximation tolerance.
tol = 1E-4;

% Load
load = -150000;

% Define step, percentage for the given load.
% 0<step<1
dF = 0.03;

% Define material properties.
E = 1E7;
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

% Initialise the internal force vector, the tangent stiffness matrix and
% th displacement vector with zeros.
qt = zeros(2*Np, 1);
kt = zeros(2*Np, 2*Np);
d = zeros(2*Np, 1);

% First prediction: internal force vector and tangent stiffness for zero
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
end

% Reduced tangent matrices(for the active DOF's).
qr = qt(afg);
kr = kt(afg,afg);

% Initialise an array to store monitor the vertical displacement of the
% loaded node and the load factor, lambda.
displacement = zeros(ceil(1/dF), 1);
la = zeros(ceil(1/dF), 1);

% Path following loop.
for i = 1:ceil(1/dF);
    
    la(i)=i*dF;
    resid=qr-la(i)*vload;
    r=norm(resid);
    
    % Approximate the target load of the current step.
    while r>tol;
        
        % Calculate the displacements for active DOFs (reduced matrix).
        dd = -kr\resid;
        d(afg) = d(afg)+dd;
        
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
            [q_e, k_e] = isoplnonlin(E, nu, th, xy, d(vu));
            
            % Add them to the global matrices.
            qt(vu) = qt(vu) + q_e;
            kt(vu,vu) = kt(vu,vu) + k_e;
        end
        
        % Reduced matrices.
        qr = qt(afg);
        kr = kt(afg,afg);
        
        % Recalculate residual.
        resid=qr-la(i)*vload;
        r=norm(resid);
    end
    displacement(i) = d(82);
end

% Global displacements.
d(afg) = d(afg)+dd;

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
plot([0; -displacement], -[0; la]*load, '-o')

