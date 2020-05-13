function [Ks, meanKs] = PermKsat (filename, L, D, method)
%PERMEAMETER calculates the water saturated hydraulic conductivity (Ks) of soil
% samples given the flux, sample length, and ponding depth (b). The
% function calculates Ks for both the constant and the falling head method.
%
% The function is based on a standard MSExcel file in which the user
% inserts flux information in three different periods of time, which requires 
% recording six time stamps).
%
%% MS Excel file structure.
%   CONSTANT HEAD (example for the first period). FIRST SHEET.
%       Tank0 = Initial water level in tank [cm].
%       Sampleh0 = 	Initial water level in the sample [cm].
%       Vol0 = 	Water level in the burette at time t0 [cm3].
%       t0 = initial time of the first period. Fraction of a day [d].
%
%       Tank1 = Final water level in tank [cm].	
%       Sampleh1 = 	Final water level in the sample [cm].
%       Vol1 = 	Water level in the burette at time t1 [cm3].
%       t1 = Final time of the first period [d].
%
%   FALLING HEAD (example for the first period). SECOND SHEET.
%       Tank0 = Initial water level in tank [cm].
%       Sampleh0 = 	Initial water level in the sample [cm].
%       t0 = initial time of the first period. Fraction of a day [d].
%
%       Tank1 = Final water level in tank [cm].	
%       Sampleh1 = 	Final water level in the sample [cm].
%       t1 = Final time of the first period [d].
%
%% Inputs
%       filename = name of the MS Excel file with the data. Follow the
%       standard file "DataPermeameter.xlsx"
%       L = Sample or ring length [cm].
%       D = ring diameter [cm].
%       method = 'constant' for constant head method.
%                'falling' for falling head method.
%
%% Outputs
%       Ks = saturated hydraulic conductivity [cm/h].
%       meanKs = mean Ks for each sample.
%

% Andres Patrignani. 29-Aug-2013 18:07:18.

%% Start code

switch method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'constant'
%% Get data
[num] = xlsread(filename);

%% Extract data into variables
% Ring Area
A = pi * (D/2).^2;
% Volume
Vol = num(:,[10 18 26]) - num(:,[6 14 22]);

% Tank mean height(h)
tank_h(:,:,1) = num(:,[4 12 20]);
tank_h(:,:,2) = num(:,[8 16 24]);
tankAVG_h = mean(tank_h,3);

% Sample mean height (h)
sample_h(:,:,1) = num(:,[5 13 21]);
sample_h(:,:,2) = num(:,[9 17 25]);
sampleAVG_h = mean(sample_h,3);

% Calculate mean ponding depth (b)
h = abs(tankAVG_h-sampleAVG_h);

% Time [d]
timemat = num(:,7:4:end);
[~, col] = size(timemat);
startvec = 1:2:col-1;
endvec = 2:2:col;
t = timemat(:,endvec)-timemat(:,startvec);

%% Calculate Ks
      Ks = (Vol*L) ./ (h.*A.*t); % [cm/d] = cm3*cm / cm*cm2*d
      meanKs = mean(Ks,2);
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'falling'   
%% Get data
[num] = xlsread('\Soil Physics Theory\rawDataPermeameter.xlsx','falling_head');

%% Extract data into variables
% Tank mean
tank_h(:,:,1) = num(:,[4 10 16]);
tank_h(:,:,2) = num(:,[7 13 19]);
tankAVG_h = mean(tank_h,3);

% Sample mean
sample_h(:,:,1) = num(:,[5 11 17]);
sample_h(:,:,2) = num(:,[8 14 20]);
sampleAVG_h = mean(sample_h,3);

% Calculate mean ponding depth (b)
 h1 = abs(tank_h(:,:,1)-sample_h(:,:,1)); % initial ponding depth for each fo the three periods.
 h2 = abs(tank_h(:,:,2)-sample_h(:,:,2)); % final ponding depth for each fo the three periods.

% Time
timemat = num(:,6:3:end);

x = 0.0864; % evaporation factor in cm/d.

%% Calculate Ks
      [~, col] = size(timemat); % take dimension of time matrix.
      startvec = 1:2:col-1; % vector with column numbers corresponding to inital time stamps.
      endvec = 2:2:col; % vector with column numbers corresponding to final time stamps.
      t = timemat(:,endvec)-timemat(:,startvec); % calculate total time for each period. 
      % By default Matlab calculates time as a fraction of a day. 
      Ks = L./t .* log ( h1./h2 ) + x*L ./ sqrt(h1.*h2) ; % Calculate Ks. Section 6.2 page 10.
      meanKs = mean(Ks,2); % Mean Ks from the three Ks values.
end

end