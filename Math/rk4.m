function [t_out,y_out] = rk4(fun,t_rng,yo,n)
%%RK4 approximates the solution of first order ODE IVPs using 4th order Runge-Kutta's method.
%
% Inputs
%       fun   = function handle of f(t,y(t)).
%       t_rng = vector containing time limits [to tf].
%       yo    = initial condition.
%       n     = number of points in which the interval to-tf will be
%               breaken into.
%
% Outputs
%       t_out = vector of n equally spaced values from to to tf.
%       y_out = vector of n values where y_out(1) = yo, and y_out(k)
%               approximates y(t) at t_out(k) for k going from 2 to n.
%
% Example 1     
%       fun = @(t,y)(y-1).^2.*(t-1).^2;
%       t_rng = [0 1];
%       yo = 0;
%       n = 20;
%
%       [t_out,y_out] = rk4(fun,t_rng,yo,n);
%
% Example 2
%       fun = @(t,y) t.*y+y+t-cos(y);
%       t_rng = [0 1];
%       yo = 1;
%       n = 20;
%
%       [t_out,y_out] = rk4(fun,t_rng,yo,n);
%
%
% Create by Andres Patrignani - 17-May-2014 17:35:33.
%
% References
% Function developed following the lecture on Runge_Kutta and Dormance Prince
% methods taught by Douglas Harder. University of Waterloo, Ontario, Canada.
% Link: http://youtu.be/Vs2usBZjUO8

% Divide t_rng into n equally spaced values.
to = t_rng(1);
tf = t_rng(2);
t_out = linspace(to,tf,n)';

% Set y_out vector and insert initial condition.
y_out = [yo; nan(n-1,1)];

% Define increment h.
h = (tf -to)/(n-1);

% Evaluate ODE (fun) at points in t_out.
for k = 1:n-1
    K1 = fun(t_out(k) , y_out(k)); % K1 is the slope at point (tk,yk) (information from first point (a) of the interval a-b
    K2 = fun(t_out(k) + h/2 , y_out(k) + h/2*K1); % See that we are using t_out(k) + h/2 [which is the same as t_out(k+h/2)] instead of t_out(k+1) as in Heun's method.
    K3 = fun(t_out(k) + h/2 , y_out(k) + h/2*K2);
    K4 = fun(t_out(k) + h , y_out(k) + h*K3);
    y_out(k+1) = y_out(k) + h*(K1 /6 + K2 /3 + K3 /3 + K4 /6); % Approximation of y at time tk+1.
    % Weighted average of K, coefficients are 1/6, 1/3, 1/3, and 1/6.
end

% Use ode45 Matlab's solver to compare our solution using Runge-Kutta 4th order method.
% ode45 uses Dormand-Prince method.
[t45,y45] = ode45(fun,t_rng,yo);

% Plot solution
figure
plot(t_out,y_out,'or');hold on;
plot(t45,y45,'-k');
legend('4th order Runge-Kutta approximation','Actual solution')
xlabel('t','FontSize',14);
ylabel('y','FontSize',14);
legend boxoff