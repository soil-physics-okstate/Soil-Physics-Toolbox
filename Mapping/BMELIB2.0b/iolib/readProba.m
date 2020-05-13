function [cs,isST,softpdftype,nl,limi,probdens,filetitle]=readProba(datafile);

% readProba                 - Reads probabilistic soft data from files (Jan 1, 2001)
%
% Reads probabilitic soft data from a proba data file.
%
% SYNTAX :
%
% [cs,isST,softpdftype,nl,limi,probdens,filetitle]=readProba(datafile);
%
% INPUT :
%
% datafile     string          filename of the data file to read. The data file
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
% OUTPUT :
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
%
%
% EXAMPLE :
% 
% Following is a valid BME proba data file:
%
%BeginningOfFile--------------------------------------------------
%BME Probabilistic data
%7    NOTE: The actual number of columns varies from line to line
%s1
%s2
%code for the variable (always equal to one)
%Type of soft pdf (always equal to 1, corresponding to histogram)
%number of limit values, nl
%limits of intervals (nl values)
%probability density (nl-1 values)
%  1  0.9   1  1  4    0.1  0.3  0.7  1.1     1.0  1.5  0.5 
%0.1  0.2   1  1  2    0.1  0.3               5.0 
%EndOfFile---------------------------------------------------------


fid=fopen(datafile,'r');
if fid==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',datafile));
end;
filetitle=fgetl(fid);
str=fgetl(fid);
nc=sscanf(str,'%d',1);
for ic=1:nc,
   strcell{ic}=fgetl(fid);
end;
if strcmp(strcell{nc-5},'time') | strcmp(strcell{nc-5},'t')
  isST=1;
else
  isST=0;
end
d=nc-5-isST;
idxcoord=1:d+isST;
idxcode=d+isST+1;
idxspt=d+isST+2;
idxnl=d+isST+3;
linenumber=nc+2;
is=0;
tempval=[];
while ~feof(fid)
  str=fgetl(fid);
  tempval=sscanf(str,'%g');
  tempval=tempval(:)';
  linenumber=linenumber+1;
  is=is+1;  
  if length(tempval)<idxnl, 
    fclose('all');
    error(sprintf('line %d of file %s should have at least %d values',linenumber,datafile,idxnl)); 
  end;
  cs{1}(is,1:d+isST)=tempval(idxcoord);
  cs{2}(is,1)=tempval(idxcode);
  if is==1
     softpdftype=tempval(idxspt);
  elseif softpdftype~=tempval(idxspt),
     fclose('all');
     error(sprintf('On line %d of file %s, the softpdftype should be the same as in previous lines',linenumber,datafile)); 
  end;
  nl(is,1)=tempval(idxnl);  
  switch softpdftype
    case {1}, nlimi=nl(is);nprobdens=nl(is)-1;
    case {2}, nlimi=nl(is);nprobdens=nl(is);
    case {3}, nlimi=3;nprobdens=nl(is)-1;
    case {4}, nlimi=3;nprobdens=nl(is);
    otherwise,
      fclose('all');
      error(sprintf('line %d of file %s softpdftype must equal to 1, 2, 3 or 4',linenumber,datafile)); 
  end;
  nval=d+isST+3+nlimi+nprobdens;  
  if length(tempval)<nval, 
    fclose('all');
    error(sprintf('line %d of file %s should have at least %d values',linenumber,datafile,nval)); 
  end;
  limi(is,1:nlimi)=tempval(d+isST+3+1:d+isST+3+nlimi);
  if softpdftype==3 | softpdftype==4
    limivec=limi(is,1):limi(is,2):limi(is,3);
    if nl(is)~=length(limivec) 
      fclose('all');
      error(sprintf('on line %d of file %s nl=%d, but limi would requires that nl=%d',...
        linenumber,datafile,nl(is),length(limivec))); 
    end;
  end
  probdens(is,1:nprobdens)=tempval(d+isST+3+nlimi+1:d+isST+3+nlimi+nprobdens);
end;
if is==0,
   cs{1}=zeros(0,d+isST);
   cs{2}=zeros(0,1);
   softpdftype=1;
   nl=zeros(0,1);
   limi=zeros(0,1);
   probdens=zeros(0,1);
end;
st=fclose(fid);
if st==-1,
  warning(sprintf('could not close %s properly',datafile));
end;

