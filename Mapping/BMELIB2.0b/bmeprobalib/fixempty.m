function [chf,zhf,csf,nlf,limif,probdensf]=fixempty(ck,ch,zh,cs,nl,limi,probdens);

% fixempty                  - Fix empty variables (Jan 1, 2001))
%
% Check if variables are empty and fix then if needed
%
% SYNTAX :
%
% [chf,zhf,csf,nlf,limif,probdensf]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
%
% INPUT :
% 
% ck         nk by d       matrix of coordinates for the estimations point,
%                          where d is the dimension of the space
% ch         nh by d or 0  matrix of coordinates for the hard data
% zh         nh by 1 or 0  vector of the values for the hard data
% cs         ns by d or 0  optional matrix of coordinates for the soft data
% nl         ns by 1 or 0  optional vector nl for soft data
% limi       ns by l or 0  optional matrix limi for soft data
% probdens   ns by p or 0  optional matrix probdens for soft data
%
% INPUT :
% 
% chf        If ch is non empty then chf=ch, else chf is a 0 by d vector
% zhf        If zh is non empty then zhf=zh, else zhf is a 0 by 1 vector
% csf        If cs is non empty then csf=cs, else csf is a 0 by d vector
% nlf        If nl is non empty then nlf=nl, else nlf is a 0 by 1 vector
% limif      If limi is non empty then limif=limi, else limif is a 0 by 1 vector
% probdensf  If probdens is non empty then probdensf=probdens, 
%            else probdensf if a 0 by 1 vector
%
% NOTE :
%
% It is possible to specify an additional index for ck,
% ch and cs, taking integer values from 1 to nv. This index
% specifies which variable is known at each coordinate.
% In that case, ck, ch, cs are cell arrays, where the first
% cell is the matrix of coordinates and the second cell is
% the vector of index values. 

if nargin<4, 
  if iscell(ck)
    cs={[],[]};
  else
    cs=[]; 
  end;
end;
if nargin<5, nl=[]; end;
if nargin<6, limi=[]; end;
if nargin<7, probdens=[]; end;


chf=ch;
csf=cs;
if iscell(ck)
  [nk d]=size(ck{1});
  if d==0, error('ck cannot be empty'); end;
  if (~iscell(chf) & isempty(chf)) | (iscell(chf) & isempty(chf{1}))
		chf{1}=zeros(0,d);
    chf{2}=zeros(0,1);
	end
  if (~iscell(csf) & isempty(csf)) | (iscell(csf) & isempty(csf{1}))
    csf{1}=zeros(0,d);
    csf{2}=zeros(0,1);
  end
else
  [nk d]=size(ck);
  if d==0, error('ck cannot be empty'); end;
  if isempty(chf)
    chf=zeros(0,d);
  end
  if isempty(csf)
    csf=zeros(0,d);
  end
end

if isempty(zh), zhf=zeros(0,1); else zhf=zh; end
if isempty(nl), nlf=zeros(0,1); else nlf=nl; end
if isempty(limi), limif=zeros(0,1); else limif=limi; end
if isempty(probdens), probdensf=zeros(0,1); else probdensf=probdens; end
