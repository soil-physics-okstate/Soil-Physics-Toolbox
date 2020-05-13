function [x]=checkmemory(pred)
%%% check_memory function - checks whether OS is Unix/Windows, then checks
% if memory is available for kriging routine (krig.m)
%
% Input:
%       pred: array of points to be predicted in kriging
%
% Output:
%       x: logical; 1 = enough memory , 0 = not enough memory
%
% Author: Matthew Haffner 12 Aug 2015
%

if ispc == 1 % The check memory function is only availble on Windows
    memavail = memory; % Check available memory
    if ((length(pred))^2)*8 > memavail.MaxPossibleArrayBytes % Check amount of memory to be used vs. amount available
        x = 0; % Not enough memory available
    else
        x = 1; % Enough memory available
    end
else
    x = 0 % If the system is not a pc, return 0
end