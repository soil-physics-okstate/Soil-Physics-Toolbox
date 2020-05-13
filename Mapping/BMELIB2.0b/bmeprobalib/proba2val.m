function [pdf]=proba2val(z,softpdftype,nl,limi,probdens)

% proba2val                 - Calculates the value of probabilistic soft pdf (Jan 1, 2001)
%
% Calculates the pdf value corresponding to different to values of the
% variable characterized by soft probabilistic data
%
% SYNTAX :
% 
% [pdf]=proba2val(z,softpdftype,nl,limi,probdens)
%
% INPUT :
% 
% z           1 by n      vector of values of the variable characterized by soft 
%                         probabilistic data.
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         probabilitic soft data.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% nl          ns by 1     vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
% limi        ns by l     matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
% probdens    ns by p     matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
%
% OUTPUT :
%
% pdf         1 by n      vectors with the value of the pdf corresponding to the 
%                         values z of the variable


if length(nl)>1 error('length(nl) must be equal 1'); end;
[zsort,idxsort]=sort(z);
pdfsort=[];
[a,b]=proba2interval(softpdftype,nl,limi);
switch softpdftype
case 1,
  idxr=find( zsort<limi(1) );
  il=1;
  while il<nl & length(pdfsort)<length(zsort)
    idxr=find( limi(il)<=zsort & zsort<limi(il+1) );
    pdfsort(idxr)=probdens(il);
    il=il+1;
  end
  pdfsort(length(pdfsort)+1:length(zsort))=0;
case 2,
  % idxr=find( zsort<limi(1) );
  % il=1;
  % while il<nl & length(pdfsort)<length(zsort)
  %  idxr=find( limi(il)<=zsort & zsort<limi(il+1) );
  %   if length(idxr)>0
  %     slope=(probdens(il+1)-probdens(il))/(limi(il+1)-limi(il));
  %     intercept=probdens(il)-slope*limi(il);
  %     pdfsort(idxr)=intercept+slope*zsort(idxr);
  %   end
  %   il=il+1;
  % end
  % pdfsort(length(pdfsort)+1:length(zsort))=0;
  pdfsort=interp1(limi,probdens,zsort);
  pdfsort(isnan(pdfsort))=0;
case 3,
  pdfsort=0*zsort;
  idxlu=find(limi(1)<zsort & zsort<limi(3));
  pdfsort(idxlu)=probdens(1+floor((zsort(idxlu)-limi(1))/limi(2)));
case 4,
  pdfsort=0*zsort;
  idxlu=find(limi(1)<zsort & zsort<limi(3));
  idxluz=1+floor((zsort(idxlu)-limi(1))/limi(2));
  limivec=limi(1):limi(2):limi(3);
  slope=(probdens(idxluz+1)-probdens(idxluz))./(limivec(idxluz+1)-limivec(idxluz));
  intercept=probdens(idxluz)-slope.*limivec(idxluz);
  slope=reshape(slope,size(idxluz));
  intercept=reshape(intercept,size(idxluz));
  pdfsort(idxlu)=intercept+slope.*zsort(idxlu);
end

pdf=pdfsort(idxsort);
pdf=reshape(pdf,size(z));
 