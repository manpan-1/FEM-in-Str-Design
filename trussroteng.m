function [q_e, K_e] = trussroteng(coor, disp, E_modulus, A_o)
%trussroteng function returning thei nternal load vector and stiffness
% matrix of a truss element using rotated engineering strain formulation.
%
% example:
% [q_e, K_e] = trussroteng(coor, disp, E_modulus, A_o)
%
% where:
% - [coor] is an array of the form [x_1, y_1; x_2, y_2] with the global 
% nodal coordinates of the element.
% - [disp] is the displacement vector of the form [u_1, v_1, u_2, v_2].
% - [E_modulus] is the modulus of elasticity.
% - [A_o] is the initial cross sectional area.

% Dimensions of the element (bounding box on the global csys).
x21=coor(3) - coor(1);
y21=coor(4) - coor(2);

% Deformations on the global csys (relative displacement of nodes).
u21=disp(3) - disp(1);
v21=disp(4) - disp(2);

% Initial and final length.
L_o = sqrt(x21^2 + y21^2);
L_n = sqrt((x21+u21)^2 + (y21+v21)^2);

% Calculate rotated engineering strain.
eps = (L_n - L_o) / L_o;

% Stress for the rotated engineering strain.
stress = E_modulus * eps;

% After the application of virtual work equilibrium, the following vector b
% is derived. Vector b multiplied with the virtual displacements gives
% virtual strains, ?? = b(d)*?d
b = 1/(L_o*L_n)*[-(x21+u21), -(y21+v21), (x21+u21), (y21+v21)];

% Calculation of the internal force vector.
q_e=A_o*L_o*stress*b';

% Auxiliary arrays for the formulation of the tangent stiffness matrix.
c = [-(x21+u21), -(y21+v21), (x21+u21), (y21+v21)];

A = [ 1,  0, -1,  0;
      0,  1,  0, -1;
     -1,  0,  1,  0;
      0, -1,  0,  1];

% Tangent stiffness matrix.
K_e = (A_o*E_modulus/(L_o*L_n^2) - stress*E_modulus/L_n^3) * (c'*c) + A_o*stress*A/L_n;

