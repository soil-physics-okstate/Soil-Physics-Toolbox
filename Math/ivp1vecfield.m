function [K,X,Y,fun_out] = ivp1vecfield (fun,t_rng,yo,meshlimits,vecsize,headsize,varargin)
%%IVP1VECFIELD generates the vector field and the solution of a first order IVP.
%
% Inputs
%       fun   = function handle of f(t,y(t)).
%       t_rng = vector containing time limits [to tf].
%       yo    = initial condition.
%       meshlimits = 
%       vecsize = 
%       headsize = 
%       n     = number of points in which the interval to-tf will be
%               breaken into.
%       method = IVP evaluation method. Possible methods are:
%                   * euler
%                   * heun
%                   * runge-kutta
%
%
%
% Outputs
%
%
%



% Example 1

% Define variables from user inputs.

%   fun = @(t,y) y.*(2-t).*t+t-1; % h is already the first derivative.

Xa = meshlimits(1,1);
Xdelta = meshlimits(1,2);
Xb = meshlimits(1,3);
Ya = meshlimits(2,1);
Ydelta = meshlimits(2,2);
Yb = meshlimits(2,3);

[X,Y] = meshgrid(Xa:Xdelta:Xb,Ya:Ydelta:Yb);
K = fun(X,Y); % slopes.
[row,col]=size(Y);
U = ones(row,col)*vecsize;
figure
quiver(X,Y,U,K,'MaxHeadSize', headsize);hold on;


% Plot a possible solution
to = t_rng(1);
tf = t_rng(2);
[t45,y45] = ode45(fun,[to tf],yo);
plot(t45,y45,'-r','LineWidth',2);
fun_out = [t45 y45];

%%% Vector analysis in 1D.
%% Example 1
syms v % Symbolic math.
f = @(v) 9.8*v - 0.098*v.^2; % Define function (falling object example)
[X,Y] = meshgrid(0:0.25:5,0:3:80); % time matrix (X) and velocity matrix (Y). 
g = gradient(f,v); % First derivative.
G = subs(g, v, Y); % Symbolic evaluation of f with respect to v, at points in Y.
[row,col]=size(Y);
U = ones(row,col).*2;
figure
quiver(X,Y,U,G,0.25,'MaxHeadSize', 0.1)
xlabel('Time');
ylabel('Velocity (m s^{-1})');

%% Example 2
syms y % Symbolic math.
h = (y.^2-y-2).*(1-y).^2; % h is already the first derivative.
[X,Y] = meshgrid(0:0.25:5,-2:0.1:3.5);
H = subs(h, y, Y);
[row,col]=size(Y);
U = ones(row,col).*20;
figure
quiver(X,Y,U,H,'MaxHeadSize', 2);hold on;

%% Example 3
syms y x % Symbolic math.
dy = y-x; % h is already the first derivative.
[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
J = subs(dy, {x,y}, {X,Y});
[row,col]=size(Y);
U = ones(row,col).*2;
figure
quiver(X,Y,U,J,'MaxHeadSize', 2);hold on;









