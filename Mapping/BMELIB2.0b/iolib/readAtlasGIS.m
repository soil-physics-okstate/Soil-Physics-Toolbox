function [pVal]=readAtlasGIS(filename);

% readAtlasGIS              - reads GIS files in Atlas GIS format (Jan 1, 2001)
%
% Attemps to reads the line objects of GIS files in Atlas GIS format. 
% (these files usually have a filename with a .bna extension).
% This program might not work for complicated data files.
%
% SYNTAX :
%
% pVal=readAtlasGIS(filename,maxline);
%
% INPUT :
%
% filename     string      name of the file to read (usually with a .bna extension)
% maxline      scalar      optional input specifying the maximum number of lines
%                          to read. Default is 10000
%
% OUTPUT :
%
% pVal        cell array   cell array of lines. Each line is a np by 2 matrix of np 
%                          2D coords, where np changes from one line to the next.
%
% EXAMPLE :
%
% Assuming you have a (small) file named test.bna in Atlas GIS format, we may be 
% able to extract the line objects and plot them in MATLAB as follow:
% 
% pVal=readAtlasGIS('test.bna');
% figure;hold on;
% for i=1:length(pVal)
%   plot(pVal{i}(:,1),pVal{i}(:,2))
% end


if nargin<2, maxlines=10000; end;

fid=fopen(filename);
i=0;
maxlines=1000;
s=fgetl(fid);                            % Get the header for the first line
while feof(fid)==0 & i<maxlines & ~isempty(s)
  comas=findstr(s,',');
  i=i+1;
  np=str2num(s(comas(end)+1:end));
  for ip=1:np
    pVal{i}(ip,1:2)=fscanf(fid,'%g, %g',2);
  end;
  fgetl(fid);
  s=fgetl(fid);                        % Get the header for the next line
end;

if i==maxlines
  warning('maximum number of lines was reached');
end

fclose(fid);
