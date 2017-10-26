function [x_c] = NewtRaph(f, y_t, x_c, tol)
%NewtRaph Newton-Raphson iteration.
% Approximates the input equation y=f(x) to a target value of
% 'y' and returns the approximated value of 'x'.
%
% example:
% x = NewtRaph([symfunc], [y_target], [x_init], [tolerance])
%
% where:
% - [symfun] is the equation to be approximated. It must be given as a
% symbolic function object (symfun) with a single symbolic variable, f(x).
% - [y_target] is the value of f(x) for which x is approximated
% - [x_init] is the initial x value to start the N-R iteration
% - [tolerance] is the approximation tolerance for the delta y

% Calculate the derivative of the given equation
eq_prime = diff(f);

% Fetch the value y of the function for the given initial x and store it as
% the current y (y_c);
y_c = double(f(x_c));

% Calculate the difference of the target y to the current y
delta_y = y_t - y_c;

% Loop through approximating steps until the differencies are smaller than
% the requested tolerance.
while tol<abs(delta_y);
    
    % Predict the next value of x based on the tangent of the current step
    % and make the result to be the new current value
    x_c = x_c + delta_y / double(eq_prime(x_c));
    
    % Calculate the value of the equation for the new current x
    y_c = double(f(x_c));
    
    % Calculate the difference of the new current y to the target y
    delta_y = y_t - y_c;
end;
