% Tutorial about implementing and understanding the Euler's method to
% approximate IVP.
dY = @(t,y) (y-1).^2 .* (t-1).^2; % ODE equation.
Y = @(t) (t.^3 - 3*t.^2 + 3*t) ./ (t.^3 - 3*t.^2 + 3*t + 3); % Known solution to the ODE.

%% Plotting actual solution using y.
t = 0:0.001:0.2;
plot(t,Y(t));hold on;

%% Plotting approximate solutions using Euler's method
% Approximate the y function at y(0.1) using information at y(0)=0.
y(1) = 0; % y0.
h = 0.1; 
yapprox = y(1) + dY(0,0)*h;
yactual = Y(y(1)+h);
error = abs(yapprox - yactual);
scatter(y(1)+h,yapprox,'or','filled');hold on;

display(['The approximate value of the function at t=',num2str(y(1)+h),...
         ' is ',num2str(yapprox)]);
display(['The actual value of the function at t=',num2str(y(1)+h),...
         ' is ',num2str(yactual)]);
display(['Using a step of h=',num2str(h),' the error is ',num2str(error(1))]);

%% Approximate the y function at y(0.05) using information at y(0)=0...
h(2) = 0.05;
yapprox = y(1) + dY(0,0)*h(2);
yactual = Y(y(1)+h(2));
error = abs(yapprox-yactual);
scatter(y(1)+h(2),yapprox,'og','filled');hold on;

display(' ');
display(['The approximate value of the function at t=',num2str(y(1)+h(2)),...
         ' is ',num2str(yapprox)]);
display(['The actual value of the function at t=',num2str(y(1)+h(2)),...
         ' is ',num2str(yactual)]);
display(['Using a step of h=',num2str(h(2)),' now the error is ',num2str(error),...
         ', which is about 4 times smaller (Remember that the error term is related',...
          ' to h square, which is the second term of the Taylor expansion that we do NOT know']);
%...and since we still want to estimate the value at t=0.1 (or to+h, which
%is 0+0.1=0.1) we need to move from the point (0.05,0.05) to the point
%(0.1,?).
%%
yapprox(2) = yapprox(1) + dY(y(1)+h(2),yapprox(1))*h(2);
error(2) = abs(yapprox(2) - Y(y(1)+h(1)));
scatter(y(1)+h(1),yapprox(2),'og','filled');hold on;

display(' ');
display(['The approximate value of the function at t=',num2str(y(1)+h(1)),...
         ' is ',num2str(yapprox(2))]);
display(['The actual value of the function at t=',num2str(y(1)+h(1)),...
         ' is ',num2str(yactual(1))]);
display(['Using a step of h=',num2str(h(2)),' the error is ',num2str(error(2))]);


