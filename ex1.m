% Exercise 1: Path following with Newton-Raphson iterations

% Give input values of the equation
z = 10;
E = 210000;
A = 1000;
L = 5000;

% Define the F(w) equation as a symbolic function
syms F(w);
F(w) = (E*A/L^3)*(z^2*w + (3/2)*z*w^2 + (1/2)*w^3);

% Define the range and the step deltaF
F_steps = 0:10:60;

% Define the tolerance for the N-R iterations approximation
tol = 0.01;

% Calculate the displacements, w, for the given force steps, F
% A function takes over the task of returning x values for a given array of
% y=f(x) values. see help steps_NewtRaph
w = steps_NewtRaph(F, F_steps, tol);

% Plot :
% a) returned displacements, w, against the corresponding F steps,
% b) F(w),w pairs for 0<w<35 with step 1
plot(w, F_steps, '-o', [0:35], F([0:35]))
legend('N-R iteration steps', 'F(w) for 0<w<35 with step 1')
xlabel('Displacement, w')
ylabel('Force, F(w)')