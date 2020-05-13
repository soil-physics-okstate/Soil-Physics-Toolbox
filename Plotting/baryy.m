function [] = baryy(x,y1,y2,width,offset,colors)
%% Generates bar plot of two series that share the same X axis.
% x = 1996:1:2014;
% y1 = rand(length(x),1);
% y2 = rand(length(x),1);
% width = 0.3;
% offset = width/2;
% colors = {'b','g'};
% baryy(x,y1,y2,width,offset,colors)

% Created by Briana Sallee and Andres Patrignani - 1 Apr 2015

if nargin == 3
    width = 0.25;
    offset = width/2;
    colors = {'b','g'};
elseif nargin == 4
    offset = width/2;
    colors = {'b','g'};
elseif nargin == 5
    colors = {'b','g'};
elseif nargin > 6
    error('Too many input arguments. Maximum number is six (6)')    
end
    
figure
plotyy(x-offset,y1,x+offset,y2,@(x,y) bar(x,y,width,colors{1}),@(x,y) bar(x,y,width,colors{2}));
