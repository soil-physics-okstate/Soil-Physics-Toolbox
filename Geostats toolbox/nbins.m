function [D,N,edges] = nbins(X)
%%NBINS Determines the optimal number of regular spaced histogram bins.
%
% Input
%       X = random variable. Column vector.
%
% Outputs
%       D = optimal number of bins
%
%       N = Number of elements in each bin
%
%       edges = boundaries of each bin
%
% Adapted by Andres Patrignani (23-Jun-2015) from Birge and Rozenholc. 2006. How many 
% bins should be put in a regular histogram. ESAIM: PS ESAIM: Probability
% and Statistics. February 2006, Vol. 10, p. 24–45.
% Defined "edges" in the loop by Jingnuo Dong 7-14-2015

% Define anonymous functions for Akaikes Information Criterion (AIC)
LogL = @(D,N,n) N*log(D*N'/n); % Eq 1.1 in Birge and Rozenholc. 2006, p.26.
pen = @(D) D-1 + log(D)^2.5; % Eq 1.2 in Birge and Rozenholc. 2006, p.26.

% Re-scale variable in the range [0,1]
newmin = 0;
newmax = 1;
Xn = rescalevar(X,newmin,newmax); % Set X in range [0,1]

% Generate a range of tentative number of bins from 1 to Dmax
% The code calculates AIC for each tentative number of bins.
n = length(X);
Dmax = round(n/log(n));
Dvec = 1:Dmax;

N = cell(Dmax,1);
edges = cell(Dmax,1);
Lpen = nan(Dmax,1);
for i=1:Dmax    
    
    %% Equation 1.1 in Birge and Rozenholc. 2006.
    I = floor(Dvec(i)*Xn)+1; % Find to which bin each values belongs to;
    [C,~,ic] = unique(I); % Get unique bins (C is incremental in each iteration 
                          % since we are testing different bins size)
    for j=1:length(C)
        count(j) = sum(ic==C(j)); % Count how many values of X are within each bin.  
    end
    
    % Merge last two bins.
    count(end-1) = count(end-1)+count(end); 
    N{i} = count(1:end-1);  % N is Nj in Birge and Rozenholc. 2006.
    edges{i} = C;
    L = LogL(Dvec(i),N{i},n); % LogLikelihood
    
    %% Eq 1.2 in Birge and Rozenholc. 2006.
    penalty = pen(Dvec(i)); % Penalty as a function of the number of bins.
    Lpen(i,1) = L-penalty; % Subtract penalty
end

% Find optimal number of bins.
[~,I] = max(Lpen); 

% Set outputs
D = Dvec(I);
edges = edges{I};
edges = rescalevar(edges,min(X),max(X));
N = N{I};