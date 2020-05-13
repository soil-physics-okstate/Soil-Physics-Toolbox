function []=writeProba(cs,isST,softpdftype,nl,limi,probdens,filetitle,datafile);

% writeProba                - Writes probabilistic soft data to files (Jan 1, 2001)
%
% Writes probabilitic soft data to a data file.
%
% SYNTAX :
%
% writeProba(cs,isST,softpdftype,nl,limi,probdens,filetitle,datafile);
%
% INPUT :
%
% cs          cell array       an array of two cells, where
%                              cs{1} is ns by d+isST matrix of coordinates, and
%                              cs{2} is ns by 1 vector of variable number
% isST        scalar           1 if last column of cs is time, 0 otherwise
% softpdftype scalar           integer for the type of soft pdf, representing the 
%                              probabilitic soft data. These soft pdf types are as follow:
%                              1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                              (see also probasyntax)
% nl          ns by 1          vector of number of interval limits. nl(i) is 
%                              the number of interval limits used to define the soft pdf 
%                              for soft data point i. (see also probasyntax)
% limi        ns by l          matrix of interval limits, where l is equal to
%                              either max(nl) or 3 depending of the softpdftype.
%                              limi(i,:) are the limits of intervals for the i-th 
%                              soft data. (see also probasyntax)
% probdens   ns by p           matrix of probability density values, where p is 
%                              equal to either max(nl)-1 or max(nl), depending on the 
%                              softpdftype. probdens(i,:) are the values of the probability 
%                              density corresponding to the intervals for the i-th soft data 
%                              defined in limi(i,:). (see also probasyntax)
% filetitle  string            title of the file
% datafile     string          filename of the data file to write. The data file
%                              has a title in the first line, the number of column headers
%                              nc in the second line, followed by nc lines with headers,
%                              followed by ns data lines. Each data line has 
%                              the coordinate of the soft data point, 
%                              a code specifying the variable number, 
%                              a code specifying the type of soft pdf (softpdftype),
%                              the number of interval limits for the variable (nl)
%                              values specifying limits of intervals for the variable (limi)
%                              and the values of the pdf at these interval limits (probdens)
%                              See probasyntax for explanation of soft probabilistic data
%
% NOTE :
%
% See also readProba

ns=size(cs{1},1);
d=size(cs{1},2)-isST;
nv=max(cs{2});

fid=fopen(datafile,'w');
fprintf(fid,'BME Probabilistic data\r',filetitle);
fprintf(fid,'%d    NOTE: The actual number of columns varies from line to line\r',d+isST+5);
switch d,   
  case 1, fprintf(fid,'s1\r');
  case 2, fprintf(fid,'s1\rs2\r');
  case 3, fprintf(fid,'s1\rs2\rs3\r');     
end;
if isST
  fprintf(fid,'time\r');
end;  
if nv>1, fprintf(fid,'code for the variable (between 1 and %d)\r',nv);
else fprintf(fid,'code for the variable (equal to one)\r'); end;
fprintf(fid,'Type of soft pdf (equal to %d, corresponding to ', softpdftype);
switch softpdftype
  case 1, fprintf(fid,'histogram)\r');
  case 2, fprintf(fid,'linear)\r');
  case 3, fprintf(fid,'histogram on regular grid)\r');
  case 4, fprintf(fid,'linear on regular grid)\r');
end;
fprintf(fid,'number of limit values, nl\r');
switch softpdftype
  case {1,2}, fprintf(fid,'limits of intervals (nl values)\r');
  case {3,4}, fprintf(fid,'limits of intervals (minval,step,maxval)\r');
end;
switch softpdftype
  case {1,3}, fprintf(fid,'probability density (nl-1 values)\r');
  case {2,4}, fprintf(fid,'probability density (nl values)\r');
end;
for i=1:ns
  str=sprintf('%15g ',cs{1}(i,:)); 
  str=[str sprintf('%3d ',cs{2}(i))];
  str=[str sprintf('%2d ',softpdftype)];
  str=[str sprintf('%4d ',nl(i))];
  switch softpdftype
    case {1,2}, str=[str sprintf('%9g ',limi(i,1:nl(i)))];
    case {3,4}, str=[str sprintf('%9g ',limi(i,:))];
  end;
  str=[str '          '];  
  switch softpdftype
    case {1,3}, str=[str sprintf('%g ',probdens(i,1:nl(i)-1))];
    case {2,4}, str=[str sprintf('%g ',probdens(i,1:nl(i)))];
  end;
  fprintf(fid,'%s\r',str);
end;
st=fclose(fid);
if st==-1,
  warning('could not close Proba data file properly');
end;

