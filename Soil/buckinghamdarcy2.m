 function [Jw,Kh,h,z] = buckinghamdarcy2 (Khmodel,modelpar,L,deltaz,h,initguess,varargin)
%%buckinghamdarcy solves the liquid water flux under unsaturated steady-state 
% flow conditions using Buckingham-Darcy equation.
%
% Inputs
%       Khmodel:  Hydraulic conductivity model (i.e. Mualem-van Genucthen, 1980)
%
%       modelpar: Model parameters can be a structure or character array.
%           *structure must have the following fields: modelpar.alpha
%                                                     modelpar.n
%                                                     modelpar.ks
%           *character array: Soil texture class (i.e. silt loam).
%
%       L:  vector containing the soil depths of the top and bottom layer.
%
%       deltaz: depth of soil layers.
%
%       h: vector describing the matric potential values at the boundaries 
%
%       initguess: initial flux (Jw) guess.
%
%       iterations: optional argument that restricts the number of
%                   iterations used to approximate the solution.
%
% Example for a silt loam soil:
%                                                      
%    alpha = 0.0051; % [1/cm].
%    n = 1.66; % [unitless]
%    Ks = 0.76; % [cm/h]
%    La = 0; % [cm]
%    ha = -1500; % [kPa]
%    Lb = 25; % [cm].
%    hb = 0; % [kPa] (water table)
%    deltaz = 0.05; % [cm].
%    Jwguess = -0.5; % [cm/h].
%    niter = 5; % number of iterations.
%    
%    modelpar.alpha = alpha;
%    modelpar.n = n;
%    modelpar.ksat = Ks; 
%
%    [Jw,Kh,h,z] = buckinghamdarcy ('vg',modelpar,[La Lb],deltaz,[ha hb],Jwguess,niter);


%% Read user inputs
if nargin == 6
    I = 10;
elseif nargin == 7
    I = varargin{1};
end

% Select hyd. cond. model
switch Khmodel
    case 'vg'
        if ~ischar(modelpar)
            alpha = modelpar.alpha; % [1/cm].
            n = modelpar.n; % [unitless].
            Ksat = modelpar.ksat; % [cm/d].
            
        elseif ischar(modelpar) % if modelpar is char find vG parameters using the RosettaParameter function.
            [~,~,~,alpha,n,~,Ksat,] = RosettaParameters (modelpar);
        else
            error('Please select a valid input argument')
        end
        
        % Define Mualem-van Genuchten, 1980 hydraulic conductivity model.
        m = 1-1/n;
        K = @(h) min(Ksat*((1-(-alpha*h).^(n-1).*(1+(-alpha*h).^n).^-m).^2) ./ (1+(-alpha*h).^n).^(m/2),Ksat);
end

%% Vertical soil column.
La = L(1); % Upper soil limit (e.g. soil surface) [cm]
Lb = L(2); % Lower soil limit [cm]
ha = h(1); % Upper boundary condition (matric potential in kPa at z = La). 
hb = h(2); % Lower boundary condition (matric potential in kPa at z = Lb).
Jw = initguess; % Initial flux guess.
z = (La:deltaz:Lb)'; % Soil depth layers.
ndepths = length(z);
h = ones(1,ndepths); % Pre-allocate matric potential array.
h(1) = ha; % Set upper soil layer boundary condition (h at z = La).

%% Iterate to find solution (Shooting method)
delta = nan(I,1); % Pre-allocate array.

for i = 1:I % Iterate to approximate Jw.
    for k = 1:ndepths-1 % Iterate to find matric potentials (h2) that satisfies the finite difference equations below.
        BuckDarcy = @(h2) Jw(i) + K( (h(k)+h2)/2 ) * ((h2-h(k))/deltaz - 1); % Set anonymous fnction in each loop to store variable values 
        h(k+1) = fzero(BuckDarcy,h(k)); % Find roots.
        
        % Error in case fzero does not converge with the initial guess.
        if isnan(h(k+1))
           error('Please choose a different initial guess');
        end
    end
    
    if i==1 % Generate a second guess by adding an arbitrary constant to the user's initial guess.
        delta(i,1) = h(end)-hb; % Calculate error from target at the lower boundary.
        Jw(i+1) = Jw(i) - 1; % Set a new flux guess (This could be replaced by the Jwmax as explained by Nofziger)
        
    else % Once we have two guesses, then we can use the Newton-Raphson method to approximate the correct flux (Jw).
         % http://en.wikipedia.org/wiki/Newton's_method
         delta(i,1) = h(end)-hb; % Evaluate difference between h obtained with current Jw and user defined lower boundary condition.
       
        % In this case the Newton-Raphson method is applied to the matric potential error (relative to the desired boundary value) as a function of Jw. 
        Jw(i+1) = interp1([delta(i-1,1) delta(i,1)],... % error between h and desired boundary value at time t and t-1.
                          [Jw(i-1) Jw(i)],...           % Flux at time t and time t-1
                          0,...                         % Tangent line x-intercept.
                          'linear','extrap');           % Extrapolate to find an improved guess of Jw. 
    end
       
end

%% Plotting commands
Kh = K(h);    

% Matric potential as a function of depth.
subplot(1,2,1),plot(h,-z);hold all;
ylabel('Distance from soil surface [cm]','FontSize',14);
xlabel('\psi_m  [kPa]','FontSize',14);
xlim([-2000 100])
 
% Hydraulic conductivity as a function of depth.
subplot(1,2,2),plot(Kh,-z);hold all;
ylabel('Distance from soil surface [cm]','FontSize',14);
xlabel('K(h)  [cm h^{-1}]','FontSize',14);
xlim([0 Ksat]);
display(['Flux (Jw) is ',num2str(Jw(i)),' cm/h'])
