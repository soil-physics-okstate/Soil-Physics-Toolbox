function [x0s,s,sv,id,b]=cokri(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%%COCKRI performs point or block cokriging in D dimensions (any integer) of
%%P variables(any integer) with combination of R basic models (any integer)
%
% Inputs
%       x: The n by (p+d) input data matrix.
%       x0: The m by d matrix of coordinates of points to estimate.
%       model: Each row of this matrix describes a different elementary 
%              structure. The first column is a code for the model type, the d 
%              following columns give the ranges along the different coordinates,
%              and the subsequent columns give rotation angles (a maximum of
%              three). The codes for the current models are:
%
%           1: nugget effect
%           2: exponential model
%           3: Gaussian model
%           4: Spherical model
%           5: Linear model
%       
%       Note: a linear model is specified by arbitrary ranges and a sill
%       such that sill/range gives the desired slope in the direction
%       considered.
%
%       c: The (rp by p) coefficient matrix of the coregionalization model.
%          Position (i,j) in each submatrix of size p by p give the sill of
%          the variable i and variable j.
%
%       itype:  Code to indicate which type of cokriging is to be perfomed:
%
%               1: Simple cokriging
%               2: Ordinary cokriging with non bias condition.
%               3: Ordinary cokriging with p nonbias condition.
%               4: Universal cokriging with drift of order 1 (linear)
%               5: Universal cokriging with drift of order 2 (quadratic)
%               99: cokriging is not performed, only sv (variances) are computed
%
%       block: Vector (1 by d), giving the size of the block to estimate;
%              any values when point cokriging is required
%
%       nd: Vector (1 by d), giving the discretization grid for block cokriging;
%           put every element equal to 1 for point cokriging.
%
%       ival: Code for cross-validation.
%
%             0: no cross-validation
%             1: cross-validation is performed by removing one variable at
%                a time at a given location.
%             2: cross-validation is performed by removing all variables at
%                a given location.
%
%       nk: Number of rearest neighboors in x matrix to use in the
%           cokriging (this includes locations with missing values even if all
%           variables are missing).
%
%       rad: Search radius for neighbors.
%
%       ntok: Points in x0 will be kriged by groups of ntok grid points.
%             When ntok>1, the search will find the nk nearest samples
%             within distance rad from the current ntok grid points
%             centroid.
%
% Ouputs
%       For the usual application, only x0s and s are required and the
%       other matrices may be ignored.
%
%       x0s: m by (d+p) matrix of the m points (blocks) to estimate by the
%            d coordinates and p cokriged estimates.
%
%       s: m by (d+p) matrix of the m points (blocks) to estimate by the d
%          coordinates and the p cokriging variances.
%
%       sv: 1 by p vector of variances of points (blocks) in the universe.
%
%       id: ((nk by p) + nc) by (ntok by p) matrix with lambda weights and
%           Lagrange multipliers of the last cokriging system solved.
%
% Coded by Andres Patrignani (minor changes were made) on 13-Jun-2015 from 
% code published by Denis Marcotte, Cokriging with matlab, Computers & 
% Geosciences, Volume 17, Issue 9, 1991, Pages 1265-1280.
%
% Function called 'means.m' in the orignal paper was not included since the
% built-in mean funciton in Matlab can do the same task.



casesen off;

% Definition of constants
[m,d] = size(x0);

% Check for cross-validation
if ival >= 1
    ntok = 1;
    x0 = x(:,1:d);
    nd = ones(1,d);
    [m,d] = size(x0);
end

[rp,p] = size(c);
[n,t] = size(x);
nk = min(nk,n);
ntok = min(ntok,m);
idp = [1:p]';
ng = prod(nd);

%--- Add d and p instead of echo

% Compute point (ng=1) or block (ng>1) variance
for i = 1:d
    n1 = prod(nd(1:i-1));
    nr = prod(nd(i+1:d));
    t = [0.5*(1/nd(i) - 1):1/nd(i):0.5*(1-1/nd(i))]';
    t2 = [t2,kron(ones(n1,1),kron(t,ones(nr,1)))]; 
end
grid = t2.*(ones(ng,1)*block);
t = [grid,zeros(ng,p)];

% For block cokriging a double grid is created by shifting slightly the
% original grid to avoid the zero distance effect (Journel and Huijbregts,
% p.96).

if ng>1
    grid = grid + ones(ng,1)*block/(ng*1e6);
end

%---l should be b according to cokri2!
[x0s,s,id,l,k0] = cokri2(t,grid,[],model,c,sv,99,avg,ng);

% sv contains the variance of points or blocks in the universe

for i = 1:p
    sv = [sv,mean(mean(k0(i:p:ng*p,i:p:ng*p),1),2)]; % Original line: sv = [sv,means(means(k0(i:p:ng*p,i:p:ng*p))')];
end

% Start cokriging

for i = 1:ntok:m
    nnx = min(m-i+1,ntok);
    disp(['kriging points # ',num2str(i),' to ',num2str(i+nnx-1)]);
    
    % Sort x samples in increasing distance relatively to centroid of
    % 'ntok' points to krige
    
    
    centx0 = ones(n,1)*mean(x0(i:i+nnx-1,:),1); % Original line: centx0 = ones(n,1)*means(x0(i:i+nnx-1,:));
    tx = [x(:,1:d) - centx0].*[x(:,1:d)-centx0]*ones(d,1);
    [tx,j] = sort(x);
    
    % Keep samples insede radius; create an identifier of each samples and
    % variable (id)
    
    t = [];
    id = [];
    ii = 1;
    tx = [tx;nan];
    while ii<nk & tx(ii)<rad^2
        t = [t;x(j(ii),:)];
        id = [id;[ones(p,1)*j(ii),idp]];
        ii = ii+i;
    end
    
    t2 = x0(i:i+nnx-1,:);
    
    % If block cokriging discretize the block
    
    t2 = kron(t2,ones(ng,1)) - kron(ones(nnx,1),grid);
    
    % Check for cross-validation
    
    if ival>=1
        est = zeros(1,p);
        sest = zeros(1,p);
        
        % Each variables is cokriged in its turn
        
        if ival == 1
            np = 1;
        else
            np = p;
        end
        
        for ip = 1:np:p
            % Because of the sort, the closest sample is the sample to
            % cross-validate and its value is in row 1 of t; a temporary
            % vector keeps the original values before perfomrming the
            % cokriging.
            
            vtemp = t(1,d+ip:d+ip+np-1);
            t(1,d+ip:d+ip+np-1) = ones(1,np)*nan;
            [x0ss,ss] = cokri2(t,t2,id,model,c,sv,itype,avg,ng);
            est(ip:ip+np-1) = x0ss(ip:ip+np-1);
            sest(ip:ip+np-1) = ss(ip:ip+np-1);
            t(1,d+ip:d+ip+np-1) = vtemp;
        end
        
        x0s = [x0s;[t2,est]];
        s = [s;[t2,sest]];
    else
        [x0ss,ss,id,l] = cokri2(t,t2,id,model,c,sv,itype,avg,ng);
        x0s = [x0s;[x0(i:i+nnx-1,:),x0ss]];
        s = [s;[x0(i:i+nnx-1,:),ss]];
    end
end


