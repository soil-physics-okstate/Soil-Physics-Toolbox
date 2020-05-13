function [csC,nlC,limiC,probdensC,idxFailed] = probacombinedupli(cs,softpdftype,nl,limi,probdens,method)

% probacombinedupli               - combine duplicate soft data of probabilistic type
%
%  Each set of duplicated soft data (with same space/time coordinates) 
%  are combined into a soft datum.  Using multiplicative combination,
%  the combined soft pdf fS(x) at a point with n duplicated soft pdfs fS1(x)
%  to fSn(x) is fS(x)=A*fS1(x)*...*fSn(x), where A is the normalization 
%  constant, whereas using the additive method it is
%  fS(x)=A*(fS1(x)+...+fSn(x))
%  The returned coordinates csC do not have any duplicated soft data point.
%
%  SYNTAX :
% 
%  [csC,nlC,limiC,probdensC,idxFailed] = probacombinedupli(cs,softpdftype,nl,limi,probdens,method)
% 
%  INPUT :
%
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
%  method      char array defining the method to combine the soft data
%              'mult'     multiplicative combination of the soft pdf
%              'add'      additive combination of the soft pdf
%              'auto'     attemps a multiplicative combination, and if it yields a 
%                         zero probdens, then performs an additive combination of pdfs
%              default is 'mult'
%
%  OUTPUT :
%
%   csC          nsC by d  matrix of the coordinates of nsC unique soft data points 
%                          in a space of dimension d. Note that nsC <= ns 
%   nlC          nsC by 1  vector of the number of interval limits.
%   limiC        nsC by l  matrix of interval limits
%   probdensC    nsC by p  matrix of probability density values
%   idxFailed    1 by n    cell array of the duplicate points for which the multiplicative 
%                          soft pdf was zero (failed multiplicative combination),
%                          so that an additive combination  was used instead. 
%                          Note that idxFailed is empty when using method 'add'
%                          
%                          
% NOTE:
%
% This function was only implemented for the case of softpdftype 
% equal to 2.  Contact the developer of BMElib if you wish to extend
% this function so it works for other softpdftype.

nlmax=100;  % Maximum value authorized for nl (see FORTRAN functions in mvnlib)

if nargin<5, error('probacombinedupli needs 5 input arguments'); end
if nargin<6, method='mult'; end
if ~isnumeric(cs), error('cs must be numeric'); end;
[ns,d]=size(cs);
softpdftypeCheckArgs(softpdftype,nl,limi,probdens);
if size(nl,1)~=ns error('cs and nl must have the same number of rows'); end

if ~isdupli(cs)
  csC=cs;
  nlC=nl;
  limiC=limi;
  probdensC=probdens;
  idxFailed={};
else 
  [iu,ir]=finddupli(cs);
  [csC,cs2,nlC,limiC,probdensC,nl2,limi2,probdens2]=probasplit(cs,softpdftype,nl,...
    limi,probdens,iu);
  nsC=length(iu);
  idxFailed={};
  for k=1:length(ir)
    nsC=nsC+1;
    csC(nsC,1:d)=cs(ir{k}(1),:);
    switch softpdftype
      case 2
        limiK=[];
        for j=1:length(ir{k})
          limiJ=limi(ir{k}(j),1:nl(ir{k}(j)));
          deltaLimiJ=min([sqrt(eps) diff(limiJ)/100]);
          probdensJ1=probdens(ir{k}(j),1);
          probdensJend=probdens(ir{k}(j),nl(ir{k}(j)));
          if probdensJ1>0, limiJ=[limiJ(1)-deltaLimiJ limiJ]; end
          if probdensJend>0, limiJ=[limiJ limiJ(end)+deltaLimiJ]; end
          limiK=[limiK limiJ];
        end
        limiK=sort(unique(limiK));
        nlK=length(limiK);
        probdensKJ=ones(size(limiK));
        for j=1:length(ir{k})
          [dummy,dummy,nlJ,limiJ,probdensJ]=probasplit(cs,softpdftype,...
            nl,limi,probdens,ir{k}(j));
          probdensKJ(j,:)=proba2val(limiK,softpdftype,nlJ,limiJ,probdensJ);
        end
        switch method
          case 'mult'
            probdensK=prod(probdensKJ,1);
            if sum(probdensK)==0
              warning('One of the combined soft pdf is zero, switching method to ''auto''');
              method='auto';
              probdensK=sum(probdensKJ,1);
              idxFailed{length(idxFailed)+1}=ir{k};
            end
          case 'add'
            probdensK=sum(probdensKJ,1);
          case 'auto'
            probdensK=prod(probdensKJ,1);
            if sum(probdensK)==0
              probdensK=sum(probdensKJ,1);
              idxFailed{length(idxFailed)+1}=ir{k};
           end
          otherwise
            error('method must be equal to ''mult'', ''add'' or ''auto''');
        end
        idx=(probdensK>0 | [0 probdensK(1:end-1)] | [probdensK(2:end) 0]);
        limiK=limiK(idx);
        nlK=length(limiK);
        probdensK=probdensK(idx);
        [probdensK]=proba2probdens(softpdftype,nlK,limiK,probdensK);
        if nlK>nlmax
          [zquantiles]=quantile(limiK',probdensK',[0.001:0.998/(nlmax-1):0.999]);
          probdensK=proba2val(zquantiles,softpdftype,nlK,limiK,probdensK);
          limiK=zquantiles;
          nlK=length(limiK);
        end
        [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,...
          softpdftype,nlK,limiK,probdensK);
      otherwise
        error('You must use softpdftype=2, contact BMElib developers if you want this function to work for other softpdftype');
    end   
  end
end

