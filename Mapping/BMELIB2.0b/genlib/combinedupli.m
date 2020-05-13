function [chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,method,cs,softpdftype,nl,limi,probdens)

% combinedupli               - combine duplicate hard and soft data
%
%  Each set of duplicated hard and soft data (with same space/time
%  coordinates) are combined into a hard or soft datum.
%  First the hard data is combined by replacing each set of duplicate hard 
%  data by either a unique hard datum (such as the average, minimum, maximum 
%  or quantile of the set of duplicated hard values), or by a soft pdf (such
%  as a student T or the histogram of the duplicated hard data).
%  Then if there are soft data, each duplicated set containing one hard datum 
%  and one or more soft data is replaced by the hard datum, and each set of 
%  duplicated soft data is combined using probacombinedupli.m
%  The returned coordinates chC and csC do not have any duplicated data point.
%
%  SYNTAX :
% 
%  [chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,method,cs,softpdftype,nl,limi,probdens)
% 
%  INPUT :
%
%  ch          nh by d    matrix of the coordinates of nh hard data points 
%                         in a space of dimension d
%  zh          n by 1     vector of zh hard data values
%  method      char array or numeric value, defining the method to combine
%                         each set of duplicated hard data as follow
%              'ave'      - average value
%              'min'      - minimum value 
%              'max'      - maximum value 
%              q          - where 0<q<1, q-quantile of each set of duplicated hard data
%              'studentT' - creates a student T soft pdf for each set of 
%                           duplicated hard data using probaStudentT.m
%              'hist'     - creates a histogram soft pdf for each set of 
%                           duplicated hard data using histscaled.m
%              default is 'ave', which is the same as using avedupli.m
%  cs          ns by d    matrix of the coordinates of ns soft data points 
%                         in a space of dimension d
%                         default is []
%  softpdftype scalar     indicates the type of soft pdf representing the  
%                         probabilitic soft data at the coordinates specified in cs.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
%                         default is 2
%   nl          ns by 1   vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
%                         default is []
%   limi        ns by l   matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
%                         default is []
%   probdens    ns by p   matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
%                         default is []
%
%  OUTPUT :
%
%   chC          nhC by d  matrix of the coordinates of nhC unique hard data points 
%                          in a space of dimension d. Note that nhC <= nh
%   zhC          nhC by 1  vector of hard data values at chC
%   csC          nsC by d  matrix of the coordinates of nsC unique soft data points 
%                          in a space of dimension d. Note that (nhC+nsC) <= (nC+nC) 
%   softpdftypeC scalar    indicates the type of soft pdf of combined soft data
%   nlC          nsC by 1  vector of the number of interval limits.
%   limiC        nsC by l  matrix of interval limits
%   probdensC    nsC by p  matrix of probability density values
%
% NOTE:
%
% softpdftype must be equal to 2 when method is equal to 'studentT' or 'hist'
%
% [chC,zhC]=combinedupli(ch,zh,'ave') returns the same output variables as
% [pu,za] = avedupli(ch,zh)

if nargin<3, method='ave'; end
if nargin<4, cs=[]; end
if nargin<5, softpdftype=2; end
if nargin<6, nl=[]; end
if nargin<7, limi=[]; end
if nargin<8, probdens=[]; end

nbins=20;  % number of bins used when method is equal to 'hist'

if ~isnumeric(ch) | ~isnumeric(cs), error('ch and cs must be numeric');  end;
[nh,d]=size(ch);
if size(zh,1)~=nh, error('zh must have the same number of rows as ch'); end
if nh~=0 & size(zh,2)~=1, error('zh must have one column'); end
if ~isnumeric(method) & ~ischar(method), error('method must be numeric or character'); end
if isnumeric(method)
  if length(method)~=1, error('when method is numeric, it must have a length of 1');
  elseif ~(0<method & method<1), error('when method is numeric, it must be such that 0<method<1');
  else 
    P=method; 
    method='quantile';
  end
end
[ns]=size(cs,1);
if ns~=0 & nh~=0 & size(cs,2)~=d, error('cs must have the same number of columns as ch'); end
if isempty(softpdftype), softpdftype==2; end;
softpdftypeCheckArgs(softpdftype,nl,limi,probdens);
if size(limi,1)~=ns, error('limi must have the same number of rows as cs'); end
if size(probdens,1)~=ns, error('probdens must have the same number of rows as cs'); end
switch method
  case {'ave','min','max','quantile'}
    methodtype='hard';
  case 'studentT'
    methodtype='soft';
    if softpdftype~=2, 
      error('When method is equal to ''studentT'' then softpdftype must be equal to 2'); 
    end
  case 'hist'
    methodtype='soft';
    if softpdftype~=2, 
      error('When method is equal to ''hist'' then softpdftype must be equal to 2'); 
    end
  otherwise
    error('Bad value for method');
end

csC=[];
softpdftypeC=softpdftype;
nlC=[];
limiC=[];
probdensC=[];

%%%%%% Combine the hard data
if nh==0 | ~isdupli(ch)
  chC=ch;
  zhC=zh;
else
  [iu,ir]=finddupli(ch);
  chC=ch(iu,:);
  zhC=zh(iu);
  nhC=length(iu);
  nsC=0;  
  for k=1:length(ir)
    switch methodtype
      case 'hard'
        nhC=nhC+1;
        chC(nhC,:)=ch(ir{k}(1),:);
        switch method
          case 'ave', zhC(nhC,1)=mean(zh(ir{k}));
          case 'min', zhC(nhC,1)=min(zh(ir{k}));
          case 'max', zhC(nhC,1)=max(zh(ir{k}));
          case 'quantile', zhC(nhC,1)=quantest(zh(ir{k}),P);
        end        
      case 'soft'
        nsC=nsC+1;
        csC(nsC,:)=ch(ir{k}(1),:);
        switch method
          case 'studentT', 
            zave(nsC,1)=mean(zh(ir{k}));
            zvar(nsC,1)=var(zh(ir{k}));
            nobvs(nsC,1)=length(zh(ir{k}));
          case 'hist',
            [n,x]=histscaled(zh(ir{k}),nbins);
            dx=diff(x(1:2));
            limi1=[x(1)-dx x' x(end)+dx];
            probdens1=[0 n' 0];
            nl1=length(limi1);
            [probdens1]=proba2probdens(softpdftype,nl1,limi1,probdens1);
            if nsC==1
              nlC=nl1;
              limiC=limi1;
              probdensC=probdens1;
            else
              [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,...
                softpdftype,nl1,limi1,probdens1);
            end
        end
    end
  end
  switch method
    case 'studentT'
      [softpdftype,nlC,limiC,probdensC]=probaStudentT(zave,zvar,nobvs);
  end
end

if isempty(nl) 
  return; 
else
  csC=[csC;cs];
  [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,...
    softpdftype,nl,limi,probdens);
  if ~isdupli([chC;csC]), return; end;
end

%%% Remove each set of duplicate soft data collocated with a hard datum
nhC=size(chC,1);
[iu,ir]=finddupli([chC;csC]);
idxRemove=[];
for k=1:length(ir)
  if min(ir{k})<=nhC
    idxRemove=[idxRemove;ir{k}(ir{k}>nhC)-nhC];
  end
end

[csRemove,csC,nlRemove,limiRemove,probdensRemove,nlC,limiC,probdensC]=probasplit(csC,softpdftype,nlC,...
  limiC,probdensC,idxRemove);

[csC,nlC,limiC,probdensC,idxFailed] = probacombinedupli(csC,softpdftype,nlC,limiC,probdensC,'auto');

if length(idxFailed)>0
  str1='Some duplicated soft data points resulted in zero soft pdf when using multiplicative combination.';
  str2='Use probacombinedupli on the soft data alone to see if there is something wrong with the soft data.';
  str3='If using probacombinedupli does not reveal any failed multiplicative combination then the problem';
  str4='might be with combining hard AND soft duplicates.  For these problematic duplicate data points we';
  str5=' had to use an additive combination of soft data (see probacombinedupli for more explanations).';
  disp(sprintf('warning in combinedupli: \n%s\n%s\n%s\n%s\n%s',str1,str2,str3,str4,str5));
end
