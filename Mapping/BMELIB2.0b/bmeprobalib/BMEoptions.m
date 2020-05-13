function options=BMEoptions

% BMEoptions                - Default parameters used in the BME functions (Jan 1, 2001)
%
% The funtions using these parameters are:
%   BMEprobaMode, BMEprobaMoments, BMEprobaPdf, BMEprobaCI
%
% options = BMEoptions returns the default options.
%
% The Parameters are
%            options(1), 1 displays the point number. (Default: 0)
%            options(2) and (14) are those used by fminbnd for finding
%                       the mode of the posterior distribution. The default
%                       values are the same as for fminbnd
%            options(3) is the maximum number of evaluation of the integral.
%                       The default value is 50000.
%            options(4) is the relative error on the estimation of the integral.
%                       The default value is equal to 1e-3.
%            options(6) is used when z=[] to determine how fine of a grid of z 
%                       values to construct. Increase options(6) to construct a finer grid. 
%                       The default value is 25
%            options(7) is used when z=[] to determine the extend of the grid of z
%                       values to construct. Decrease options(7) to widen the range of
%                       the grid of z values. The default value is 0.001
%            options(8) number of moments to calculate (1, 2 or 3)
%                       =1 to calculate the mean of the BME posterior pdf,
%                       =2 for the mean and estimation error variance,
%                       =3 for the mean, estimation error variance and coef of skewness.
%                       The default value is 2
%            options(10) this options is reserved for the return argument of the fminbnd function
%            options(14) see options(2)
%            options(20:29) confidence probability for which to compute the confidence interv
%                       The CI are calculated only for options(20:29) which are >=0.01 and <=0.99
%                       The default value is options(20)=0.68 and options(21:29)<0
%            
%
%  BMEprobaMode uses options 1, 2, 3, 4, 14
%  BMEprobaMoments uses options 1, 3, 4, 8
%  BMEprobaPdf uses options 1, 2, 3, 4, 6, 7, 14
%  BMEprobaCI uses options 1, 2, 3, 4, 6, 7, 14, 20:29
% 
%            The other values of options are not used
options(1)=0;
options(2)=1e-4;
options(3)=50000;
options(4)=1e-3;
options(6)=25;
options(7)=0.001;
options(8)=2;
options(14)=100;
options(20)=0.68;
options(21)=0;
options(22)=0;
options(23)=0;
options(24)=0;
options(25)=0;
options(26)=0;
options(27)=0;
options(28)=0;
options(29)=0;

    
