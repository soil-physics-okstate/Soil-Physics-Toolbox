function TC = devries(theta_l,T,theta_i,rho_b,P,dg,lambda_s,q0,gs)
%% deVries
% Calculates the thermal conductivity of porous media. 
%
%% Syntax
%
% |TC = devries(theta_l,T,theta_i,rho_b,P,dg,lambda_s,q0,gs)|
%
%% Inputs
%       Parameter   Variable                        Units
%       theta_l     liquid water content            m3 m-3
%       T           temperature                     K
%       theta_i     ice content                     m3 m-3
%       rho_b       bulk density                    kg m-3
%       P           atmospheric pressure            Pa
%       dg          geometric mean grain diameter   microns
%       lambda_s    thermal conductivity of solids  W m-1 K-1
%       q0          liquid recirculation exponent   n/a
%       gs          shape factor for soil solids    n/a
%
%% See also
% <deVriesSaturationExample.html |deVriesSaturationExample|>
%
%% References
%
% de Vries, D.A. 1963. Physics of Plant Environment by W.R. van Wijk. Chapter 7.
%
% Campbell et al. 1994. Soil Sci. 158:307-313.
%
% Ochsner, T.E., Horton, R. and Ren, T.H. 2001. A new perspective on soil 
% thermal properties. Soil Sci Soc Am J 65: 1641-1647.
%
%% Updates
%
% Created by Tyson E. Ochsner, November 30, 2006

%%
%

% Step 1:  Determine the continuous medium
vs = rho_b/2650;  % volume fraction of solids (2650 kg m3 is the density of the solids)
theta_s = 1 - vs; % saturated water content
vg = 1 - vs - theta_l - theta_i; % air content

media = [0 0 0 0];
if theta_l > 0.2*theta_s
    media(1)=1; %water is the continuous media
elseif theta_i > vg
    media(2)=1; %ice is the continuous media
else
    media(3)=1; %air is the continous media
end
media = logical(media);

% Step 2:  Calculate thermal conductivities of the various constituents
lambda_w = 0.554 + 2.24*10^-3*(T-273.15)-9.87*10^-6*(T-273.15)^2; %water
lambda_a = 0.024 + 7.73*10^-5*(T-273.15)-2.6*10^-8*(T-273.15)^2; %dry air
lambda_i = 2.16 - 1.01*10^-2*(T-273.15)+5.28*10^-5*(T-273.15)^2; %ice

% The following are needed to calculate the thermal conductivity of the gas
% phase according to the procedure of Campbell et al. 1994. Soil Science
% 158:307-313.

Hv = 45144 - 48*(T-273.15); % latent heat of vaporization
h = 1; % relative humidity
q = q0*(T/303)^2; % exponent for water recirculation function
xw0 = 0.267*dg^(-0.2); % cutoff water content for liquid recirculation
fw = [1+(theta_l/xw0)^-q]^-1; % weighting function from Campbell et al. 1994
P0 = 101325; % atmospheric pressure at sea level (Pa)
rho_0 = 44.65; %standard density of dry air (mol m-3)
rho_hat = rho_0*(P/P0)*(273.15/T); %density of dry air
Dv0 = 2.12*10^-5; % standard diffusivity of water vapor in air (m2 s-1)
Dv = Dv0*(P0/P)*(T/273.15)^1.75; % diffusivity of water vapor in air
t = 1 - 375.15/T; % dimensionless temperature
p_star = P0*exp(13.3016*t - 2.042*t^2 + 0.26*t^3 + 2.69*t^4); % sat vap press (Pa)
s = 375.15*p_star*(13.3016 - 4.082*t + 0.78*t^2 +10.76*t^3)/T^2; % slope

lambda_v = Hv*h*fw*rho_hat*Dv*s/(P-h*p_star); % latent heat effect on lambda
lambda_g = lambda_a + lambda_v; % effective thermal conduvtivity of gas phase

lambdaarray = [lambda_w lambda_i lambda_g lambda_s];
lambda_0 = lambdaarray(media); % select the conductivity of the continous media

%Step 3:  Calculate the weighting factors

%Solids
ks = 2/3*[1+(lambda_s/lambda_0 - 1)*gs]^-1 + 1/3*[1+(lambda_s/lambda_0 - 1)*(1-2*gs)]^-1;

%Water
gw = gs; %As assumed by de Vries (1963) and Hopmans and Dane (1986)
kw = 2/3*[1+(lambda_w/lambda_0 - 1)*gw]^-1 + 1/3*[1+(lambda_w/lambda_0 - 1)*(1-2*gw)]^-1;

%Gas
%The shape factor function below is my simplification based on the work of
%de Vries (1963), Hopmans and Dane (1986), and Campbell et al. (1994)
if theta_l > xw0
    gg = 1/3;
else
    gg = 1/3*(theta_l/xw0);
end
kg = 2/3*[1+(lambda_g/lambda_0 - 1)*gg]^-1 + 1/3*[1+(lambda_g/lambda_0 - 1)*(1-2*gg)]^-1;

%Ice
%Ice would tend to follow the same pore geometry as air, so the shape
%factor is the same
if theta_l > xw0
    gi = 1/3;
else
    gi = 1/3*(theta_l/xw0);
end
ki = 2/3*[1+(lambda_i/lambda_0 - 1)*gi]^-1 + 1/3*[1+(lambda_i/lambda_0 - 1)*(1-2*gi)]^-1;

%Step 5:  Calculate thermal conductivity

TC = (kw*theta_l*lambda_w + ki*theta_i*lambda_i + ks*vs*lambda_s + kg*vg*lambda_g)/...
    (kw*theta_l + ki*theta_i + ks*vs + kg*vg);