function [q_e, k_e, so_el, epso_el] = isoplnonlinplast(th, x, d, so_el, epso_el, do)
% ISOPLNONLIN calculates and returns the internal force vector and tangent
% stiffness matrix for a non-linear 2D planar element.
%
% Input variables
% th : shell thickness
% x  : nodal coordinates
% d  : nodal displacements
% so_el : previous steps stresses of the element (all Gauss points)
% epso_el : previous step plastic stress of th element
% do : previous step nodal displacements

% Get globals
global E nu

% Calculate diferential strain based on the previous and the current steps'
% displacement
dd = d - do;

% Assemble elasticity matrix.
% C_m = E/(1-nu^2)*[1,  nu, 0;
%     nu, 1,  0;
%     0,  0,  (1-nu)/2];

% Matrix addressing the diferential displacements to the 3 strains ?x, ?y,
% ?xy.
H_m = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];

% Initialise internal force vector and tangent stiffness matrix.
q_e = zeros(8, 1);
k_e = zeros(8,8);

% Gauss points local coordinates.
a = 1/sqrt(3);
xi = [-a -a  a a];
mu = [-a  a -a a];

% Assemble stiffness matrix.
for i = 1:4
    % Jacobian
    Dn = 1/4*[-(1-mu(i)),  (1-mu(i)), (1+mu(i)), -(1+mu(i));
              -(1-xi(i)), -(1+xi(i)), (1+xi(i)),  (1-xi(i))];
    
    J = Dn*[x(1), x(2)
            x(3), x(4)
            x(5), x(6)
            x(7), x(8)];
    % Inverse of the Jacobian.
    Ga = inv(J);
    
    % Formulate the transformation matrix to to convert local to global
    % csys expressions.
    Gama_m = [Ga, zeros(2, 2)
          zeros(2, 2), Ga];
    
    N_m = 1/4*[-(1-mu(i)), 0,        (1-mu(i)), 0,       (1+mu(i)), 0,       -(1+mu(i)), 0;
              -(1-xi(i)), 0,       -(1+xi(i)), 0,       (1+xi(i)), 0,        (1-xi(i)), 0;
               0,        -(1-mu(i)), 0,       (1-mu(i)), 0,        (1+mu(i)), 0,       -(1+mu(i));
               0,        -(1-xi(i)), 0,      -(1+xi(i)), 0,        (1+xi(i)), 0,        (1-xi(i))];
    
    % G matrix (connects the global displacements vector to the
    % differential displacements on the local level)
    G_m = Gama_m * N_m;
    
    % A matrix (represents the non-linear part of B)
    A_m = [G_m(1, :)*d, 0,           G_m(3, :)*d, 0;
           0,           G_m(2, :)*d, 0,           G_m(4, :)*d;
           G_m(2, :)*d, G_m(1, :)*d, G_m(4, :)*d, G_m(3, :)*d];
    
    % Calculate the B matrix (connecting virtual displacements to virtual
    % strains).
    B_m = (H_m + A_m)*G_m;
    
    % Calculate differential strains (from diff. displacement from previous
    % step for plasticity).
    de = B_m*dd;
    
    % Current step strains.
    e_1 = (H_m + 0.5*A_m)*G_m*d;
    
    % Previous step strain
    e_0 = (H_m + 0.5*A_m)*G_m*do;
    
    % The deifferential strain de can be instead calculated as the
    % difference of the strains, e_0, e_1
    % de e_1-e_0
    
    % Calculate the new stresses using pstress2d function.
    [S_m, epsn, C_m] = pstress2d(so_el{i},epso_el(i),de);
    
    % Stress matrix.
    %S_m = C_m*e_1;
    
    % Define D matrix (stresses for tangent stiffness matrix)
    D_mm = [S_m(1), S_m(3); S_m(3), S_m(2)];
    D_m = [D_mm, zeros(2); zeros(2), D_mm];
    
    % Internal force vector.
    q_e = q_e + (B_m'*S_m);
    
    % Tangent stiffness matrix.
    k_e = k_e + th*(B_m'*C_m*B_m + G_m'*D_m*G_m);
    
    % Replace the so and epso to the newlly calculated
    so_el{i} = S_m;
    epso_el(i) = epsn;
    
end;

k_e = k_e*th;