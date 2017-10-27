function [x] = steps_NewtRaph(equation, y, tol)
% steps_NewtRaph Series of Newton-Raphson steps.
% For a given equation y=f(x) and given values y, calculate and return the
% corresponding values of x, approximated using Newton Raphson iterations.
%
% example:
% [x_array] = steps_NewtRaph([symfunc], [y_target], [tolerance])
%
% where, 
% - [symfun] is the equation to be approximated. It must be given as a
% symbolic function object (symfun) with a single symbolic variable, f(x).
% - [y_target] is the array of f(x) values for which x is approximated
% - [tolerance] is the approximation tolerance for the delta y


% Start the first N-R iteration from an arbitary value, u=0 and approximate
% the first y value, y(1)
x(1) = NewtRaph(equation, y(1), 0, tol);

% Loop through the rest of the y values and run N-R iterations. On each
% step, the N-R iteration starts from the previous approximated value of x.
for i = 2:length(y);
    x(i) = NewtRaph(equation, y(i), x(i-1), tol);
end;


