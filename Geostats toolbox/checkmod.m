function checkmod(model,c,d,rad)
% CHECKMOD: This function generates 'ntot' points within a D-sphere of
% radius 'rad' (in reduced rotated distances) and evaluates all
% cross-variograms and variograms as described in 'model' and 'c'.
%
% Inputs
%       model: description of models (see COKRI)
%       c: matrices of sills (see COKRI)
%       d: dimension (1-D, 2-D,...)
%
% Coded by Andres Patrignani on 13-Jun-2015 from code published by Denis 
% Marcotte, Cokriging with matlab, Computers & Geosciences, Volume 17, 
% Issue 9, 1991, Pages 1265-1280.


[t,p] = size(c);
[r,t] = size(model);
ntot = 3000/max(d,p*p)

% Equations for the variograms, models should be placed in the same order
% as in COKRI2 note that here the vartiograms are computed instead of the
% covariograms as in COKRi2.

Gam = ['h~=0'; % nugget
       '1-exp(-h)'; % exponential
       '1-exp(-h.^2)'; % gaussian
       '(1.5*min(h,1)/1 - 0.5*min(h,1)/1).^3)'; % spherical
       'h']; % linear

 % Generates 'ntot' random points inside a D-sphere of radius 'rad'
 
 t = rand(ntot,d) - 0.5*ones(ntot,d)*rad*2;
 
 % Evaluate variograms for each simulated point
 
 k = zeros(ntot*p,p);
 for i = 1:r
     t = trans(t,model,i);
     h = sqrt(sum([t.*t]')');
     ji = (i-1)*p+1;
     js = i*p;
     g = eval(Gam(model(i,1),:));
     k = k + kron(g,c(ji:js,:));
 end
 
 % Check that absolute value of semivariance is less than
 % sqrt(semivariance^2)
 
 for i = 1:ntot
     ii = (i-1)*p+1;
     is = i*p;
     t = sqrt(diag(k(ii:is,:)));
     k(ii:is,:) = abs(k(ii:is,:)) - t*t';
 end
 
 k = [k>1e-10];
 
 if sum(sum(k)>0)
     warning('Warning model is not admissible')
 end