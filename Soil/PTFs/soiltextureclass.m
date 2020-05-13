function [T] = soiltextureclass (sand, clay)
%%SOILCLASS retrieves soil texture classification using sand and clay percentages.
% Soil texture classes are defined by a point given by the amount of sand,
% clay, and silt.
%
%% Inputs
%   sand = col vector (n by 1). Values must be in percentages.
%   clay = col vector (n by 1). Values must be in percentages.
%
%% Outputs
%   T = n by m character array without apostrophes. n is
%   the number of soils and m is the number of characters in the longest
%   soil texture class name of the output.
%
%% Example
%
%   [T] = soiltextureclass ([20;40;30],[20;20; 35]);
%
%% References
% The following reference provides non-overlaping soil texture classes:
% E. Benham and R.J. Ahrens, W.D. 2009. Clarification of Soil Texture Class 
% Boundaries. Nettleton National Soil Survey Center, USDA-NRCS, Lincoln, Nebraska.
%
% Created: 19-Apr-2014 14:03:16 by Andres Patrignani.

%% Check user inputs
% Check for column vector.
if ~iscolumn(sand) || ~iscolumn(clay)
    error('input arguments must be column vectors');
end

% Check for equal length
if numel(sand) ~= numel(clay)
   error('input arguments must have the same number of elements')
end

% Check for entries adding more than 100%.
if sum([sand clay],2)>100
    error('one or more entries add over 100%')
end

% Check for entries lower than zero.
if ~isempty(sand(sand<0)) || ~isempty(clay(clay<0)); % Alternative: sum(sign(sand)<0)>=1 || sum(sign(clay)<0)>=1
   error('one or more entries are negative')
end

%% Generate cell array with soil texture classes. There are 12 different 
% classes according to United States Department of Agriculture (USDA).
% An array of soil texture classes (12 by n) is created to match a logical indexing matrix.
texturalClasses = repmat({'Sand',...
                          'Loamy sand',...
                          'Sandy loam',...
                          'Loam',...
                          'Silt loam',...
                          'Silt',...
                          'Sandy clay loam',...
                          'Clay loam',...
                          'Silty clay loam',...
                          'Sandy clay',...
                          'Silty clay',...
                          'Clay'}',...
                          1,length(sand)); 
       
%% Calculate silt content.    
silt = 100-sand-clay;       

%% Define boundaries for each soil texture class.
% All entries (rows) are evaluate for each class to create a logical indexing array.
class (1,:) = silt + 1.5*clay < 15; % Sand.

class (2,:) = silt + 1.5*clay >= 15 & silt + 2*clay < 30; % Loamy sand.

class (3,:) = (clay >= 7 & clay < 20 & sand > 52 & silt + 2*clay >= 30) |...
              (clay < 7 & silt < 50 & silt + 2*clay >= 30); % Sandy loam.               

class (4,:) = clay >= 7 & clay < 27 & silt >= 28 & silt < 50 & sand <= 52; % Loam.

class (5,:) = (silt >= 50 & clay >= 12 & clay < 27) |...
              (silt >= 50 & silt < 80 & clay < 12); % Silt loam.

class (6,:) = silt >= 80 & clay < 12; % Silt.

class (7,:) = clay >= 20 & clay < 35 & silt < 28 & sand > 45; % Sandy clay loam.

class (8,:) = clay >= 27 & clay < 40 & sand > 20 & sand <= 45; % Clay loam.

class (9,:) = clay >= 27 & clay < 40 & sand <= 20; % Silty clay loam.

class (10,:) = clay >= 35 & sand > 45; % Sandy clay.

class (11,:) = clay >= 40 & silt >= 40; % Silty clay.

class (12,:) = clay >= 40 & sand <= 45 & silt < 40; % Clay.

%% Set output
T = char(texturalClasses(class)); % Set output without apostrophes. Ideal for copy and paste outside Matlab.

                               