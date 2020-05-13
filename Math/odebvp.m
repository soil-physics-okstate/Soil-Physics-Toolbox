function [x_out,u_out] = odebvp (C,g,x_boundary,u_boundary,n)
%%ODEBVP solves any second or first order linear ODE with known coefficients.
%
% Inputs
%       C = vector with constants
%       g = forcing function (anonymous or .m file)
%       x_boundary = vector with baoundary conditions of independent variable x.
%       u_boundary = vector with boundary conditions of dependent variable.
%       n =  number of points to approximate the analytical solution.
% 
% Outputs
%       x_out = vector of n arbitrary values for the x variable.
%       u_out = vector of n values evaluated at x_out.
%
% Example for solving a second order linear DE with constant coefficients ( 1u''(x) + 2u'(x) + 3u(x) = 0 ).
% 
%       C = [1 3 2];
%       g = @(x)  x.*0; % Forcing function is equal to zero (homogeneous)
%       x_boundary = [0,1];
%       u_boundary = [4 5];
%       n = 20;
%
%       [x_out,u_out] = odebvp (C,g,x_boundary,u_boundary,n)
%       
%       A slightly different solution is obtained if constant are changed:
%       
%       hold all;
%       [x_out,u_out] = odebvp ([1.5 2.2 3],g,x_boundary,u_boundary,n);
%
%  This code was done following the lectures of Douglas Harder in Boundary-value 
%  Problems and Finite-difference Equations. University of Waterloo,
%  Ontario, Canada.
%
%  Created by Andres Patrignani, 03-May-2014 19:56:34


% Generate arbitrary points of the independent variable x.
x_out = linspace(x_boundary(1), x_boundary(2), n)';
h = 0.125; % Set h to an arbitrary value. h is the delta value in order to 
% approximate the first and second derivatives according to Taylor's series 
% expansion. This could be included as a function input.

% Simplify notation resulting from the finite difference equation that
% approximates a second order ODE using centered difference (see notes for
% more detail).
dminus = 2*C(1)-h*C(2); % Value that will be in the upper diagonal of matrix D.
d = 2*h.^2*C(3) - 4*C(1); % Value that will be in the main diagonal of matrix D.
dplus = 2*C(1) + h*C(2); % Value that will be in the lower diagonal of matrix D.

% Build matrix D so that D*u = f, and then we solve for u using u = D-1*f 
Dminus = diag(dminus .* ones(n-3,1),-1); % Upper diagonal.
D = diag(d .* ones(n-2,1)); % Main diagonal.
Dplus = diag(dplus .* ones(n-3,1),1); % Lower diagonal.

D =  Dminus + D + Dplus; % Merge matrices to create matrix D.

% Build vector f using g(x).
f = 2 * h.^2 .* g(x_out(2:end-1)); % 2*h^2 is the term we used to eliminate 2h and h^2 in the denominator of the LHS.
f(1) = f(1) - dminus * u_boundary(1); % First value of vector includes the known value at the u boundary.
f(end) = f(end) - dplus * u_boundary(2); % Last value of vector f includes the other value at the boundary.

% Estimate u values.
u_out = D \ f; % Inverse matrix multiplication.
u_out = [u_boundary(1); u_out; u_boundary(2)]; % Insert u boundary values.

% Plot
%plot(x_out,u_out,'or');hold on;
plot(x_out,u_out)
xlabel('x','Fontsize',16)
ylabel('u(x)','Fontsize',16)





