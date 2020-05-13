function [val,valname,filetitle]=readGeoEAS(datafile,colnum);

% readGeoEAS                - Reads data files in GEO format (Jan 1, 2001)
%
% Reads specified columns of a data files in GEO format.
%
% SYNTAX :
%
% [val,valname,filetitle]=readGeoEAS(datafile,colnum);
%
% INPUT :
%
% datafile   string   data file in Geo EAS format. This format has a title in the
%                     the first line, the number of column nc in the second line
%                     followed by nc lines with the name of each column,
%                     followed by the data lines, each data lines having nc value
% colnum     vector   optional parameter, specifying the columns to read. 
%                     max(colnum) cannot be greater than nc. Default value is 1:nc
%
% OUTPUT :
%
% val        n by nv  matrix of value read for the nv variables, where nv=sum(colnum>0)
% valname    cells    array of nv cells with the name of the nv variables (for each colum
%                     of val)
% filetitle  string   char array with the title of the file
%
% EXAMPLE :
% 
% Following is an example of a valid file in GeoAES format:
%
%BeginningOfFile--------------------------------------------------
%BME hard data
%3
%s1
%s2
%primary variable
%10.1    11.4    1.4 
%12.2    13.6    1.6 
%16.7    19.1    1.1 
%10.9    16.9    0.9 
%EndOfFile--------------------------------------------------------


fid=fopen(datafile,'r');     % Open the data file
if fid==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',datafile));
end;
filetitle=fgetl(fid);        % Read the title of the data file
str=fgetl(fid);              % Read the number of column in data file
nc=sscanf(str,'%d',1);
if ~isnumeric(nc) | length(nc)~=1
   fclose('all');
  error(sprintf('Problem reading the number of variables on line 2 of %s',datafile));
end

if nargin<2
  colnum=1:nc;
end;
if sum(colnum<=0)>1  
  error('All the column number specified in colnum must be greater than 0');
end
if length(colnum)==0
  val=[];
  valname=[];
  return;
end;

if nc<max(colnum)
  error(sprintf('%s has a request for colomn %d, but it has only %d columns',datafile,max(colnum),nc));   
end;
for ic=1:nc,                 % Read the name of the data file columns
  str=fgetl(fid);
  datafilevalname{ic}=sscanf(str,'%s',1);
end;
for i=1:length(colnum)      % Assign the name of the specified columns to valname
  valname{i}=datafilevalname{colnum(i)};
end
linenumber=nc+2;
i=0;
tempval=[];
while ~feof(fid)
  str=fgetl(fid);            % Read each data line
  tempval=sscanf(str,'%g'); 
  linenumber=linenumber+1;
  i=i+1;
  if length(tempval)~=nc, 
    fclose('all');
    error(sprintf('line %d of file %s should have %d values',linenumber,datafile,nc)); 
  end;
  tempval=tempval(:)';
  val(i,1:length(colnum))=tempval(colnum);
end;
if i==0,
   val=zeros(0,length(colnum));
   valName=[];
end;

st=fclose(fid);
if st==-1,
  warning(sprintf('could not close %s properly',datafile));
end;

