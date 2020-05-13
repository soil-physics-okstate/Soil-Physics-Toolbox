function [xscaled,xpar] = rescalevar(x,varargin)
%%RESCALEVAR re-scales a given vector to a new desired range.
%
% Syntax
%        [xscaled,xpar] = rescalevar(x,newmin,newmax)
%   
%        [xscaled,xpar] = rescalevar(x,newminmax) first row represents the 
%                         new minimum values for each column in x, and
%                         second row must contain the new maximum values 
%                         for each column in x.
%
% Example
% M = magic(4)
% [xscaled,xpar] = rescalevar(M,0,1) % Re-scale each column to range [0,1]
% N = rescalevar(xscaled,xpar)       % Recover original values
%
%
%
% Andres Patrignani - 14 Apr 2015
% AP - Added code to recover scaled values - 8 Jul 2015

if nargin == 1 || nargin > 3
    error('input:number',['Input arguments must be two or three. You had ',num2str(nargin),' input arguments'])
    
elseif length(varargin) == 1 && size(varargin{1},1) == 2 && size(varargin{1},2) == size(x,2)
    newmin = varargin{1}(1,:);
    newmax = varargin{1}(2,:);
    newmin = repmat(newmin,size(x,1),1);
    newmax = repmat(newmax,size(x,1),1);
    
elseif length(varargin) == 2 && numel(varargin{1}) == 1 && numel(varargin{2}) == 1
    newmin = varargin{1};
    newmax = varargin{2};
    newmin = ones(size(x))*newmin;
    newmax = ones(size(x))*newmax;
end

xmax = max(x);
xmin = min(x);
xpar = [xmin;xmax];

xmax = repmat(xmax,size(x,1),1);
xmin = repmat(xmin,size(x,1),1);

xscaled = (newmax-newmin)./(xmax-xmin).*(x-xmax) + newmax; 
