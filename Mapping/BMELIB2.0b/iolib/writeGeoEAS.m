function []=writeGeoEAS(val,valname,filetitle,datafile,colnum,outofrangeval);

% writeGeoEAS               - Writes data files in Geo EAS format (Jan 1, 2001)
%
% Writes data in specified columns of a Geo EAS data files.
%
% SYNTAX :
%
% writeGeoEAS(val,valname,filetitle,datafile,colnum,outofrangeval);
%
% INPUT :
%
% val           n by nv    matrix with the values to write.
% valname       cell       1 by nv cell array with the name of the 
%                          variables (one for each colum of val)
% filetitle     string     title of the file
% datafile      string     data file in Geo EAS format. This format has a title in the
%                          the first line, the number of column nc in the second line
%                          followed by nc lines with the name of each column,
%                          followed by the data lines, each data lines having nc value
% colnum        vector     optional paramater specifying the column assigned to each variable.
%                          length(colnum) must be equal to the number of column of val.
%                          max(colnum) is the number of column nc of the datafile
%                          If columns of datafile that are not assigne to a variable
%                          are filled with outofrangeval.
%                          Default value is colnum=1:size(val,2)
% outofrangeval scalar     option value to use when writing NaNs. (Default=0)
%
% NOTE :
%
% See also readGeoEAS



if nargin<5, colnum=1:size(val,2); end
if nargin<6, outofrangeval=0; end

if sum(colnum<=0)>1  
  error('All the column number specified in colnum must be greater than 0');
end
if length(colnum)~=size(val,2)
  error('length(colnum) must be equal to the number of column of val');
end;

nc=max(colnum);
nvar=size(val,2);

fid=fopen(datafile,'w');       % Open the data file
if fid==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',datafile));
end;

fprintf(fid,'%s\r',filetitle);  % Write the title of the data file
fprintf(fid,'%d\r',nc);         % Write the number of columns
for ic=1:nc,                    % Write the name of the data file columns
  ivar=find(colnum==ic);
  if isempty(ivar)
    fprintf(fid,'%s\r','column not used');
  else   
    fprintf(fid,'%s\r',valname{ivar});
  end
end;

nval=size(val,1);               % Reassign val to valfile
valfile=outofrangeval*ones(nval,nc);
valfile(1:nval,colnum)=val;
valfile(isnan(valfile))=outofrangeval;

for i=1:nval
  str=sprintf('%12.10g ',valfile(i,:)); 
  fprintf(fid,'%s\r',str);      % Write each data line
end

st=fclose(fid);
if st==-1,
  warning(sprintf('could not close %s properly',datafile));
end;

