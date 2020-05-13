function [main,sub] = customlayout(rows,cols,vspace,hspace,axsize,margin)
%CUSTOMLAYOUT returns POSITION matrices for a figure and its subplots.  The units
%of the MAIN matrix are determined by the units of the other arguements 
%which must be consistent.  Definitions of the inputs are:
% 
% rows = number of rows of subplots
% cols = number of columns of subplots
% vspace = vertical spacing between subplots
% hspace = horizontal spacing between subplots
% axsize = [width height] of subplot axes
% margin = width of margin around the outside
%
% Tyson Ochsner
% v1 March 14, 2006
% v2 February 7, 2007

axwidth = axsize(1);    %width of axis 
axheight = axsize(2);   %height of axis 
numsubs = rows*cols;    %number of subplots

figheight = axheight*rows+(rows-1)*vspace+2*margin;
figwidth = axwidth*cols+(cols-1)*hspace+2*margin;
main = [1 1 figwidth figheight];

height = axheight/figheight;          % normalized height of axes
width = axwidth/figwidth;            % normalized width of axes
vspace = vspace/figheight;          % normalized vertical spacing
hspace = hspace/figwidth;          % normalized horizontal spacing
leftini = margin/figwidth;          %left side of first subplot
bottomini = margin/figheight;       %bottom of first subplot

n = 1;

for i = 1:rows
    for j = 1:cols
        left = leftini + (j-1)*(width + hspace);
        bottom = bottomini + (rows-i)*(height + vspace);
        sub(n,:) = [left bottom width height];
        n = n+1;
    end
end

