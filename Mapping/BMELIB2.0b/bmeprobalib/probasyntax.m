function probasyntax

% probasyntax               - Syntaxical help for soft probabilistic data (Jan 1, 2001)
%
% probasyntax provides a help file to explain the syntax of soft 
% probabilistic data in BMELIB.
%
% SYNTAX :
%
% help probasyntax;
%
% or, more simply, you may just type:
%
% probasyntax;
%
%
% CONTENT OF HELP FILE:
%
% Soft probabilistic data are defined by means of soft probability
% density functions (pdf). These soft pdfs may have different types, 
% which are explained here. 
%
% Let's consider ns soft data points with corresponding vector
% of random variable Xsoft=[X1,...,Xns]. The soft pdf fs(Xsoft) is 
% then the ns-dimensional pdf fs(Xsoft)=fs1(X1)*...*fsns(Xns). 
% Along each of the ns dimension the univariate pdf fs(Xi) is 
% defined using intervals of values for Xi. The interval limits are 
% specified using the matrix limi, and the value of fsi in these 
% intervals is specified by the matrix probdens.  Hence a soft pdf 
% is defined using four variables: pdfsofttype, nl, limi, and probdens.
%
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         probabilitic soft data.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, constant pdf value in each interval
%                         2 for Linear, linear pdf value in each interval
%                         3 for Histogram on a regular grid (intervals have same size)
%                         4 for Linear on a regular grid
% nl          ns by 1     vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. 
% limi        ns by l     matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         If softpdftype==1 or 2 then limi is a ns by max(nl) matrix,
%                           and limi(i,1:nl(i)) are the interval limits for soft data i.
%                         If softpdftype==3 or 4 then limi is a ns by 3 matrix. The 
%                           interval limits are on a regular grid, and limi(i,1:3) are
%                           the lower limit, increment, and upper limit of the interval 
%                           limits for soft data i.
% probdens    ns by p     matrix of probability density values, where p is equal 
%                         to either max(nl)-1 or max(nl), depending on the softpdftype.
%                         If softpdftype==1 or 3 then probdens is a ns by max(Nl)-1 matrix.
%                           The pdf value is constant in each interval, and 
%                           probdens(i,nl(i)-1) are the value of the pdf in each interval.
%                         If softpdftype==2 or 4 then probdens is  a ns by max(Nl) matrix.
%                           The pdf value varies linearly between interval limits, 
%                           and probdens(i,nl(i)) are the value of the pdf at each 
%                           interval limit.
%
%
% IMPORTANT NOTE: A valid pdf must be such that the integral of the pdf
%                  be equal to one.
%
% EXAMPLES:
%   All the following soft pdfs represent a uniform distribution in [0,1]
%   softpdftype=1;nl=2;limi=[0,1];  probdens=[1];
%   softpdftype=2;nl=2;limi=[0,1];  probdens=[1,1];
%   softpdftype=3;nl=2;limi=[0,1,1];probdens=[1];
%   softpdftype=4;nl=2;limi=[0,1,1];probdens=[1,1];
%
% SEE ALSO: softpdftypeCheckArgs

help probasyntax;
  
