function [pVal]=readARCe00(filename,maxline)

% readARCe00                - reads ARCINFO files in .e00 format (Jan 1, 2001)
%
% Attemps to reads the line objects of a ARCINFO file in e00 format 
% (these files usually have a filename with a .e00 extension). This program
% might not work for complicated e00 files. 
%
% SYNTAX :
%
% pVal=readARCe00(filename,maxline)
%
% INPUT :
%
% filename    string       name of the file to read (usually with a .e00 extension)
% maxline     scalar       optional input specifying the maximum number of lines
%                          to read. Default is 10000
%
% OUTPUT :
%
% pVal        cell array   cell array of lines. Each line is a np by 2 matrix of np 
%                          2D coords, where np changes from one line to the next.
%
% EXAMPLE :
%
% Assuming you have a (small) file name test.e00 in ARCINFO format, we may be 
% able to extract the line objects and plot them in MATLAB as follow:
% 
% pVal=readARCe00('test.e00');
% figure;hold on;
% for i=1:length(pVal)
%    plot(pVal{i}(:,1),pVal{i}(:,2))
% end

if nargin<2, maxlines=10000; end;

fid=fopen(filename);
if fid==-1
  error(sprintf('Could not open file %s',filename));
end
fgetl(fid);
fgetl(fid);

i=1;
nv{i}=fscanf(fid,'%d',7);
n1=nv{i}(1);
while feof(fid)==0 & i<maxlines & n1~=-1
  np(i)=nv{i}(length(nv{i}));
  fgetl(fid);
  pVal{i}=reshape(fscanf(fid,'%g %g',2*np(i)),2,np(i))';
  i=i+1;
  nv{i}=fscanf(fid,'%d',7);
  n1=nv{i}(1);
  % disp('press a key to continue'); pause; disp('continuing');
end;

if i==maxlines
  warning('maximum number of lines was reached');
end

fclose(fid);
