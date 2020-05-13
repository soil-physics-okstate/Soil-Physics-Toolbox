function fsame=fcmp(file1,file2)

% fcmp                      - Compare two files (Jan 1, 2001)
%
% Compare two files, and return 1 if the files are identical,
% return 0 is the files are different, and return -1 if one of
% file could not be open.
%
% SYNTAX :
%
% [fsame]=fcmp(file1,file2)
%
% INPUT :
%
% file1   string      name of file 1
% file2   string      name of file 2
%
% OUTPUT :
%
% fsame   scalar      integer equal to: 
%                         1 if the files are identical,
%                         0 if the files are different, and 
%                        -1 if one of file could not be open.


fid1=fopen(file1,'r');
fid2=fopen(file2,'r');
if fid1==-1 
  fsame=-1;
  return;
end;
if fid2==-1
  fsame=-1;
  return;
end;
il=0;
fsame=1;
while ~feof(fid1) & ~feof(fid2) & fsame
   il=il+1;
   l1=fgetl(fid1);
   l2=fgetl(fid2);
   if ~strcmp(l1,l2) 
      fsame=0; 
   end;
   if feof(fid1) & ~feof(fid2)
      fsame=0; 
   end;
   if ~feof(fid1) & feof(fid2)
      fsame=0; 
   end;      
end;   
fclose(fid1);
fclose(fid2);
   
giveWarning=1;  % 0-does not give warning, 1-warns which line is diff
if giveWarning & ~fsame    
  if ~strcmp(l1,l2) 
     disp(sprintf('line %d of files %s and %s is different',il,file1,file2));
     disp(sprintf('file %s:\n %s',file1,l1));
     disp(sprintf('file %s:\n %s',file2,l2));
  else   
     disp(sprintf('files %s and %s do not have the same number of lines',file1,file2));
  end;
end;
