function [P,Pind,res]=contingency(z1,z2);

% contingency            - contigency (frequency) table between two categorical variables
%                          (December 1, 2003)
%
% Compute the frequencies of the joint occurences of categories 
% for two categorical variables
%
% SYNTAX : 
%
% [P,Pind,res]=contingency(z1,z2);
%
% INPUT :
%
% z1      n by 1       vector of codes for the categories of the first
%                      categorical variable. Categories are coded as
%                      integers ranging from 1 to nc1, where nc1 is the
%                      number of possible categories.
% z2      n by 1       vector of codes for the categories of the second
%                      categorical variable. Categories are coded as
%                      integers ranging from 1 to nc2, where nc2 is the
%                      number of possible categories.
%
% OUTPUT :
%
% P       nc1 by nc2   table of observed frequencies for category i of the
%                      first variable and category j of the second variable
% Pind    nc1 by nc2   table of frequencies for category i of the first
%                      variable and category j of the second variable under
%                      the hypothesis that variables are independent (joint
%                      probabilities P(i,j) are obtained as the product of
%                      marginal probabilities.
% res     nc1 by nc2   Pearson's residuals, i.e. standardized differences between
%                      observed frequencies P(i,j) and frequencies Pind(i,j)
%                      under independence hypothesis.
%
% NOTE : If independence hypothesis holds, asymptotic distribution for the
%        res(i,j) values are zero mean unit variance Gaussian distributions.

%%% Initialize parameters

ncat1=max(z1);
ncat2=max(z2);
n=length(z1);
P=zeros(ncat1,ncat2);

%%% Compute marginal probabilities

for i=1:ncat1,
  p1(i)=sum(z1==i)/n;
end;

for i=1:ncat2,
  p2(i)=sum(z2==i)/n;
end;

%%% Compute joint probabilities

for i=1:ncat1,
  for j=1:ncat2,
    P(i,j)=sum((z1==i)&(z2==j));
  end;
end;
P=P/n;

Pind=p1'*p2;

%%% Compute Pearson's residuals

res=(P-Pind)/sqrt(Pind.*(1-Pind));

