function BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax,order);

% BMEprobaCheckArgs         - Check arguments for the BMEproba functions (Jan 1, 2001)
%
% Checks the inputs arguments for BMEproba functions such as BMEprobaMode. 
% Returns an error message if a problem is found.
%
% SYNTAX :
%
% BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
%   covmodel,covparam,nhmax,nsmax,dmax,order);
%
% INPUT : Same as the BMEproba functions (such as BMEprobaMode).
%

%%%% Test types of ck, ch and cs, and get sizes
if isnumeric(ck),
  nk=size(ck,1);           % nk is the number of estimation points
  if ~isnumeric(ch) error('If ck is a matrix then ch must be a matrix too'); end;
  nh=size(ch,1);           % nh is the number of hard data
  if ~isnumeric(cs) error('If ck is a matrix then cs must be a matrix too'); end;
  ns=size(cs,1);           % ns is the number of soft data
  dk=size(ck,2);           % dk is the dimension of ck space
  dh=size(ch,2);           % dh is the dimension of ch space
  ds=size(cs,2);           % ds is the dimension of cs space
  d=dk;
  if dh~=d, error('ck and ch must have the same number of columns'); end;
  if ds~=d, error('ck and cs must have the same number of columns'); end;
  nv=1; 
  if isdupli(ch)
    error('Duplicate hard data points.  Use finddupli to find them, and combinedupli to average them');
  elseif isdupli(cs)
    error('Duplicate soft data points.  Use finddupli to find them, and probacombinedupli to combine them');
  elseif isdupli([ch;cs])
    error('Duplicate soft and hard data points.  Use finddupli to find them, and combinedupli to combine them');
  end
elseif iscell(ck)
  if ~iscell(ch) error('If ck is a cell array then ch must be a cell array too'); end;
  if ~iscell(cs) error('If ck is a cell array then cs must be a cell array too'); end;
  if min(size(ck))~=1 | length(ck)~=2, 
     error('ck must be a 2 by 1 cell array');
  end;  
  if min(size(ch))~=1 | length(ch)~=2, 
     error('ch must be a 2 by 1 cell array');
  end;  
  if min(size(cs))~=1 | length(cs)~=2, 
     error('cs must be a 2 by 1 cell array');
  end;
  nk=size(ck{1},1);
  nh=size(ch{1},1);
  ns=size(cs{1},1);
  if nk>0 & ~isnumeric(ck{1}), error('ck{1} must be numeric'); end;
  if nh>0 & ~isnumeric(ch{1}), error('ch{1} must be numeric'); end;
  if ns>0 & ~isnumeric(cs{1}), error('cs{1} must be numeric'); end;
  if nk>0 & ~isnumeric(ck{2}), error('ck{2} must be numeric'); end;
  if nh>0 & ~isnumeric(ch{2}), error('ch{2} must be numeric'); end;
  if ns>0 & ~isnumeric(cs{2}), error('cs{2} must be numeric'); end;
  dk=size(ck{1},2);   
  dh=size(ch{1},2);     
  ds=size(cs{1},2);
  d=dk;
  if dh~=d, error('ck{1} and ch{1} must have the same number of columns'); end;
  if ds~=d, error('ck{1} and cs{1} must have the same number of columns'); end;
  if size(ck{2},1)~=nk | size(ck{2},2)~=1
     error('ck{2} must be a nk by 1 vector');
  end;
  if size(ch{2},1)~=nh | size(ch{2},2)~=1
     error('ch{2} must be a nh by 1 vector');
  end;
  if size(cs{2},1)~=ns | size(cs{2},2)~=1
     error('cs{2} must be a ns by 1 vector');
  end;
  nv=max([ck{2};ch{2};cs{2}]);
  if isdupli([ch{1} ch{2}])
    error('Duplicate hard data points.  Use finddupli to find them, and combinedupli to average them');
  elseif isdupli([cs{1} cs{2}])
    error('Duplicate soft data points.  Use finddupli to find them, and probacombinedupli to combine them');
  elseif isdupli([[ch{1} ch{2}];[cs{1} cs{2}]])
    error('Duplicate soft and hard data points.  Use finddupli to find them, and combinedupli to combine them');
  end
else
  error('ck must be either a matrix or a cell of length 2');
end;

%%%%% Check the hard data
if ~isnumeric(zh) error('zh must be a numeric vector'); end;  
if size(zh,1)~=nh error('ch and zh must have the same number of rows'); end;
if nh>0 & size(zh,2)~=1, error('zh must be a nh by 1 vector'); end;

%%%%% Check the soft data
if size(nl,1)~=ns error('cs and nl must have the same number of rows'); end;
softpdftypeCheckArgs(softpdftype,nl,limi,probdens);
if ns>0
  [probdensnorm,norm]=proba2probdens(softpdftype,nl,limi,probdens);
  if min(norm)<0.95 | max(norm)>1.05
    warning('soft probabilistic data had to be normalized to be a proper pdf');
    probdens=probdensnorm;
  end;
end;

%%%%% Check the covariance model and parameter, 
if isnumeric(ck) & ischar(covmodel)                             %%%%% case for 1 variable and one model
  if ~isnumeric(covparam) 
    warning('covmodel is a char array but covparam is not a numeric array'); 
  end; 
  if min(size(covparam))~=1 
     warning('covparam is not a vector');
  end;
elseif iscell(ck) & ischar(covmodel),                           %%%%% case for nv variables and one model
  if ~iscell(covparam)
    error('When using nv variables and 1 model, covparam must be a cell array');
  end;
  if min(size(covparam))~=1 | length(covparam)~=2
    error('When using nv variables and 1 model, covparam must be a vector of 2 cells');
  end;
  if ~isnumeric(covparam{1})
    error('When using nv variables and 1 model, covparam{1} must be a numeric matrix');       
  end;
  if size(covparam{1},1)~=nv | size(covparam{1},2)~=nv
    error('When using nv variables and 1 model, covparam{1} must be a nv by nv matrix');       
  end;
  if  ~isnumeric(covparam{2})
    error('When using nv variables and 1 model, covparam{2} must be numeric');       
  end;
  if min(size(covparam{2})~=1)
    error('When using nv variables and 1 model, covparam{2} must be a vector');       
  end;   
elseif isnumeric(ck) & iscell(covmodel),                        %%%%% case for one variable and nm models
  if min(size(covmodel))~=1,
    error('when covmodel is a cell array, it must be a vector of cells')
  end;
  nm=length(covmodel);                                          % compute the number of models 
  if ~iscell(covparam)
     error('if covmodel is a cell array then covparam must be a cell array too');
  end;
  if min(size(covparam))~=1,
    error('when covparam is a cell array , it must be a vector of cells')
  end;
  if length(covparam)~=nm,
     error('covmodel and covparam must be cell arrays with same length');
  end;
  for i=1:nm,                                                   % check elements of cell
    if ~ischar(covmodel{i})
      error('when covmodel is a cell array, each of its element must be a character string');
    end;
    if ~isnumeric(covparam{i})
       error('When using 1 variable and nm models, each cell of covparam must be a numeric vector');
    end;
  end;  
elseif iscell(ck) & iscell(covmodel),                           %%%%% case for nv variables and nm models 
  if min(size(covmodel))~=1,
    error('when covmodel is a cell array, it must be a vector of cells')
  end;
  nm=length(covmodel);                                          % compute the number of models 
  if ~iscell(covparam)
     error('if covmodel is a cell array then covparam must be a cell array too');
  end;
  if min(size(covparam))~=1,
    error('when covparam is a cell array, it must be a vector of cells')
  end;
  if length(covparam)~=nm,
     error('covmodel and covparam must be cell arrays with same length');
  end;
  for i=1:nm,                       % check cell elements 
    if ~ischar(covmodel{i})
      error('when covmodel is a cell array, each of its cell must be a character string');
    end;
    if ~iscell(covparam{i})
       error('When using nv variables and nm models, for i=1:nm, covparam{i} must be a cell array');
    end;
    if min(size(covparam{i}))~=1 | length(covparam{i})~=2
       error('When using nv variables and nm models, for i=1:nm, covparam{i} must be a vector of 2 cells');
    end;
    if ~isnumeric(covparam{i}{1})
       error('When using nv variables and nm models, for i=1:nm, covparam{i}{1} must be a numeric matrix');       
    end;
    if size(covparam{i}{1},1)~=nv | size(covparam{i}{1},2)~=nv
       error('When using nv variables and nm models, for i=1:nm, covparam{i}{1} must be a nv by nv matrix');       
    end;
    if  ~isnumeric(covparam{i}{2})
       error('When using nv variables and nm models, for i=1:nm, covparam{i}{2} must be numeric');       
    end;
    if length(covparam{i}{2})>1
      if size((covparam{i}{2}),1)~=1
        error('When using nv variables and nm models, for i=1:nm, covparam{i}{2} must be a row vector');       
      end  
    end;   
  end;    
else
  error('ck must either be numeric or a cell array, while covmodel must either be a char string or a cell array');   
end;
[isST,isSTsep,modelS,modelT]=isspacetime(covmodel);

%%%%%% Check nhmax, nsmax, and dmax
if ~isnumeric(nhmax)
   error('nhmax must be numeric');
end;
if min(size(nhmax))~=1 | length(nhmax)~=nv
   error('when using nv variables, nhmax must be a vector of length nv');
end;
if ~isnumeric(nsmax)
   error('nsmax must be numeric');
end;
if min(size(nsmax))~=1 | length(nsmax)~=nv
  error('when using nv variables, nsmax must be a vector of length nv');
end;
if sum(nsmax)>20,
  error('sum(nsmax) must not exceed 20');
end;
if ~isnumeric(dmax)
   error('dmax must be numeric');
end;
if ~isST
  if length(dmax)~=1
   error('When using a purely spatial or temporal covariance model, dmax must be a scalar');
 end;
elseif isST
  if length(dmax)~=3
    error('When using a space/time covariance models, dmax must be a vector of length 3');
  end;
end;
