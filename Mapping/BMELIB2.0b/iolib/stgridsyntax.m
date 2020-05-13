% stgridsyntax              - Syntaxical help for s/t grid format (Jan 1, 2001)
%
% Provides a help file explaining the syntax of s/t grid format versus the
% s/t vector format.
%
% SYNTAX :
%
% help stgridsyntax;
%
% or, more simply, you may just type:
%
% stgridsyntax;
%
% CONTENT OF HELP FILE:
%
% When the data are collected at a fixed set of nMS Monitoring Sites with
% coordinates cMS for nME time events, then one can use the s/t grid format to
% store the space/time data. In this case the locations of the Monitoring
% Sites are stored in the nMS by 2 vector cMS=[s1MS s2MS], the time of the
% nME Monitoring Events are stored in a 1 by nME vector, and the space/time 
% data is stored in a nMS by nME matrix Z, as shown in the following example:
%
%                tME=[100   200   300  400];
% 
% cMS=[10 10;    Z=  [1     4     7    10;
%      20 20;         2     5     8    11;
%      30 30];        3     6     9    12];
%    
% One can convert the data from s/t grid format to s/t vector format 
% where each data value is associated with a s/t location. Hence the 
% values are in a nMS*nME by 1 vector z to which is associated the  
% locations ch=[s1,s2,t], where s1, s2, and t are nMS*nME by 1 vectors 
% for the s1 spatial coordinates, s2 spatial coordinates, and time, 
% respectively.
%
% For example converting the above s/t grid data to s/t vector format
% is done by using the following command :
%
% [ch,z]=valstg2stv(Z,cMS,tME);
%
% And the result is as follow in s/t vector format
% 
%     s1    s2   time        z
%     
% ch=[10    10   100;     z=[1;
%     20    20   100;        2;
%     30    30   100;        3;
%     10    10   200;        4;
%     20    20   200;        5;
%     30    30   200;        6;
%     10    10   300;        7;
%     20    20   300;        8;
%     30    30   300;        9;
%     10    10   400;        10;
%     20    20   400;        11;
%     30    30   400];       12];
%
% Additionally the data can be converted from s/t vector to s/t
% grid format as follow
% 
% [Z]=valstv2stg(ch,z,cMS,tME);
%
% Note that the commands probastv2stg and probastg2stv do a similar
% conversion on probabilistic data.

help stgridsyntax