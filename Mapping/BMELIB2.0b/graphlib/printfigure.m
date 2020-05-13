function printfigure(figurehandle,device,filename);
%  printfigure                      - print figures to a printer device
%
%  Prints the figures to a printer device
%
%  SYNTAX :
%
%  printfigure(figurehandle,device,filename);
%
%  INPUT :
% 
%  figurehandle   vector or 'all'  vector of handles of the figures to print,
%                                  or 'all' to specify that all existing figures are to be printed
%                                  default is to print only the current figure
%  device         string or cell   a string of char specifying a printer device (see command print),
%                                  or a cell array of strings specifying several printer devices
%                                  default value is '-djpeg'
%  filename       string or cell   a string of char specifying the name of the file in which the
%                                  figure is printed, or
%                                  a cell of of same length as figurehandle containing the filenames
%                                  to use for each figure handle
%                                  by default names starting with 'figure' and followed by a number
%                                  are automatically created.
%
%  EXAMPLES :
% 
%  printfigure;    %  prints the current figure to a jpeg file, with a filename automatically generated
%  printfigure('all')  % prints all the figures open to jpeg files, with filenames automatically generated


if nargin<1
  figurehandle=get(0,'Children'); 
  if length(figurehandle)==0, return; end;
  figurehandle=figurehandle(1);
end;
if ischar(figurehandle)
  if ~strcmp(figurehandle,'all')
    error('specified figure does not exist')
  else
    figurehandle=get(0,'Children'); 
    if length(figurehandle)==0, return; end;
    figurehandle=figurehandle(end:-1:1);
  end    
end

if length(figurehandle)==0, return; end;

if ~checkfigs(figurehandle), error('Invalid figure handle.'); end;

if nargin<2, device={'-djpeg'};
elseif ischar(device), device={device};
elseif ~iscell(device)
  error('device must be either a character string or a cell of character strings');
end

filenameextension=[];
if nargin<3
  switch device{1}
  case '-djpeg', filenameextension='.jpg';
  case {'-dps','-dpsc'}, filenameextension='.ps';
  case {'-deps','-depsc'}, filenameextension='.eps';
  otherwise
    warning('You need to update printfigure.m to handle the filename extension for this printer driver');
  end  
  i0=1;
  while exist(['./figure' sprintf('%2.2d',i0) filenameextension],'file')
    i0=i0+1;
  end
  for i=1:length(figurehandle)
    filename{i}=['./figure' sprintf('%2.2d',i+i0-1) filenameextension];
  end
end
  
if ischar(filename)
  filename={filename};
end
    
if length(figurehandle)~=length(filename)
  error('figurehandle and filename must have the same length');
end;

for i=1:length(figurehandle)
  figure(figurehandle(i));
  set(gcf,'PaperUnits','inches'); 
  set(gcf,'PaperPosition',[1 1 5.5 3.4]);
  switch length(device)
  case 1, print(device{1},filename{i});
  case 2, print(device{1},device{2},filename{i});
  case 3, print(device{1},device{2},device{3},filename{i});
  case 4, print(device{1},device{2},device{3},device{4},filename{i});
  otherwise,
    warning('You need to update printfigure.m to handle more than 4 device entries');  
    print(device{1},device{2},device{3},device{4},filename{i})  
  end
end


%------------------------------------------------
function status = checkfigs(h)
status = 1;
for i=1:length(h)
  if ~ishandle(h(i)) | ~strcmp(get(h(i),'type'),'figure')
    status = 0;
    return
  end
end
