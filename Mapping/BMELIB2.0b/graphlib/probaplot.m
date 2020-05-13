function [h]=probaplot(softpdftype,nl,limi,probdens,S,idx)

% probaplot                 - Plots the probabilistic data (Jan 1, 2001) 
%
% [h]=probaplot(softpdftype,nl,limi,probdens) plots the probabilistic
%     soft data and returns the handle h. See help probasyntax for 
%     explanation of the soft data defined by softpdftype, nl, limi
%     and probdens.
% 
% [h]=probaplot(softpdftype,nl,limi,probdens,S) 
%     plots the probabilistic soft data using the linetype indicated
%     by S, where S is a character string made from one element
%     from any or all the following 3 colunms:
% 
%           y     yellow        .     point              -     solid
%           m     magenta       o     circle             :     dotted
%           c     cyan          x     x-mark             -.    dashdot 
%           r     red           +     plus               --    dashed   
%           g     green         *     star
%           b     blue          s     square
%           w     white         d     diamond
%           k     black         v     triangle (down)
%                               ^     triangle (up)
%                               <     triangle (left)
%                               >     triangle (right)
%                               p     pentagram
%                               h     hexagram
% 
% [h]=probaplot(softpdftype,nl,limi,probdens,S,idx) plots only the soft
%    pdf associated with the index idx of soft data point.
% 
% OUTPUT :
% h        vector of length ns of handles to the lines representing the 
%          soft data at each soft data points, where ns=length(nl)
%

if nargin<5, S='b-'; end;
if nargin<6, idx=1:length(nl); end;

ns=length(nl);

for i=1:length(idx)
  is=idx(i);
  switch softpdftype
  case 1,
    x=kron(limi(is,1:nl(is)),[1 1]);
    y=[0 kron(probdens(is,1:nl(is)-1),[1 1]) 0];
  case 2,
    x=[limi(is,1) limi(is,1:nl(is)) limi(is,nl(is))];
    y=[0 probdens(is,1:nl(is)) 0];
  case 3,
    x=kron(limi(is,1):limi(is,2):limi(is,3),[1 1]);
    y=[0 kron(probdens(is,1:nl(is)-1),[1 1]) 0];
  case 4,
    x=[limi(is,1) limi(is,1):limi(is,2):limi(is,3) limi(is,3)];
    y=[0 probdens(is,1:nl(is)) 0];
  end
  hold on;
  h(i)=plot(x,y,S);
end;

