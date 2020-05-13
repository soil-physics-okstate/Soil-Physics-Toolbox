function [bestCombination, minVar] = grpminvar (inputVector, grpsize)
%GRPMINVAR Function that finds the combination with the minimum variance of
%a specific input vector and for a given group size.
%
%% Inputs
% |inputVector| must be a column or row vector of the variable of
% interest. This function can handle NaN by using nan var function.
% |grpsize| must be a positive integer higher than zero, and equal or lower
% than the |length of |inputVector|. If grpsize has decimal values, grpsize
% will be rounded to the nearest integer.
%
%% Output
% |bestCombination| is a row vector of size |grpsize| by n that represents 
% the combination of elements in inputVector that have the minimum variance.

% Andres Patrignani
% 03-Aug-2013 14:43:06

%% Detect user's error

if  grpsize > length(inputVector)    % Check for size of grpsize.
    error('grpsize cannot be greater than size of inputVector')
    
elseif grpsize <=0  % Check for grpsize being negative or equal to zero.
    error ('grpsize cannot be lower or equal to zero')
    
elseif ~isnumeric(inputVector) || ~isnumeric(grpsize)   % Check for numeric input.
    error('input arguments must be numeric')
    
elseif rem(grpsize, 1) ~= 0 % Check for grpsize bein an integer.
    grpsize = round(grpsize);
   
end

%% Calculate all possible combinations (nchoosek)
combos = combntns(1:numel(inputVector),grpsize);

%% Pre-allocate memory
grpVar = zeros(size(combos,2),1);

%% Calculate variance for each combination
for i = 1:size(combos)
    grpVar(i,1) = nanvar(inputVector(combos(i,:))); % nanvar is used to handle NaN.
end

%% Find value and index of minimum variance
[minVar, minVarIdx] = min(grpVar);

%% Use minVarIdx to find the combination with minimum variance
bestCombination = combos(minVarIdx,:); 

end