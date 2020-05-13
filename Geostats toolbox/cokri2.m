function [x0s,s,id,b,k0] = cokri2(x,x0,id,model,c,sv,itype,avg,ng)
%
% Function called within COKRI. The description for input and output is
% given in COKRI. THe only new variables are 'k0' which is the right member
% matrix of the cokriging system and 'ng' which is the total number of
% points for block discretization.
%
% Define semivariogram models. New model can be added here.
%
% Coded by Andres Patrignani (minor changes were made) on 13-Jun-2015 from 
% code published by Denis Marcotte, Cokriging with matlab, Computers & 
% Geosciences, Volume 17, Issue 9, 1991, Pages 1265-1280.

Gam = ['h==0'; % nugget
       'exp(-h)'; % Exponential
       'exp(-h^2)'; % Gaussian
       '1-(1.5*min(h,1)/1-0.5*(min(h,1)/1).^3)'; % Spherical
       '1-h']; % Linear
   
% Definition of constants

[n,t] = size(x);
[rp,p] = size(c);
r = rp/p;
[m,d] = size(x0);
cx = [x(:,1:d);x0];

% Calculation of left covariance matrix K and right covariance matrix K0

k = zeros(n*p,(n+m)*p);

for i = 1:r
    
    % Calculation of matrix of reduced rotated distances H
    
    [t] = trans(cx,model,i);
    t = t*t';
    h = sqrt(-2*t+diag(t)*ones(1,n+m) + ones(n+m,1)*diag(t)');
    h = h(1:n,:);
    ji = (i-1)*p+1;
    js = i*p;
    
    % Evaluation of the current basic structure
    
    g = eval(Gam(model(i,1),:));
    k = k + kron(g,c(ji:js,:));
end
k0 = k(:,n*p+1:(n+m)*p);
k = k(:,1:n*p);

% Constraints are added according to cokriging type

if itype == 99
    
    % The system does not have to be solved
    return
end

if itype == 2
    
    % Cokriging with one non-bias condition (Isaaks and Srivastava, 1990,
    % p.410.)
    
    k = [k,ones(n*p,1); ones(1,n*p),0];
    k0 = [k0;ones(1,m*p)];
    nc =1;
elseif itype>=3
    
    % Ordinary cokriging (Myers, Math. Geol, 1982)
    
    t = kron(ones(1,n),eye(p));
    k = [k,t';t,zeros(p,p)];
    k0 = [k0;kron(ones(1,m),eye(p))];
    nc = p;
    
    if itype>=4
        
        % Universal cokriging; linear drift constraints
        
        t = kron(cx(1:n,:),ones(p,1));
        k = [k,[t;zeros(p,d)];[t',zeros(d,d+p)]];
        t = kron(cx(n+1:n+m,:)',ones(1,p));
        k0 = [k0;t];
        nc = nc+d;
    end
    
    if itype==5
        
        % Universal cokriging; quadratic drift constraints
        
        nca = d*(d+1)/2;
        cx2=[];
        for i = 1:d
            for j = i:d
                cx2 = [cx2,[cx(:,i).*cx(:,j)]];
            end
        end
        t = kron(cx2(1:n,:),ones(p,1));
        k = [k,[t;zeros(nc,nca)]; [t',zeros(nca,nc+nca)]];
        t = kron(cx2(n+1:n+m,:)', ones(1,p));
        k0 = [k0;t];
        nc = nc + nca;
    end
end

% Columns of k0 are summed up (if necessary) for block cokriging

m = m/ng;
t = [];

for i = 1:m
    for ip = 1:p
        j = ng*p*(i-1)+ip;
        
        t = [t,mean(k0(:,j:p:i*ng*p),2)']; % Original line: t = [t,means(k0(:,j:p:i*ng*p)')'];
    end
end
k0 = t;
%end NOT SURE IF THIS END IS CORRECT OR IS AN EXTRA END

t = x(:,d+1:d+p);

if itype<3
    
    % If simple cokriging or cokriging with one non bias condition, the
    % means are subtracted.
    
    t = (t-ones(n,1)*avg)';
else
    t = t';
end

% Removal of lines and columns in K and K0 corresponding to missing values

z = zeros(n*p,1);
z(:) = t;
iz = ~isnan(z);
iz2 = [iz;ones(nc,1)];
nz = sum(iz);

% If no samples left, return NaN

if nz == 0
    x0s = nan;
    s = nan;
    return
else
    k = k(iz2,iz2');
    k0 = k0(iz2,:);
    id = id(iz,:);
    
    % Solution of the cokriging system by gaussian elimination
    l = k\k0;
    
    % Calculation of cokriging estimates
    
    t2 = l(1:nz,:)'*z(iz);
    t = zeros(p,m);
    t(:) = t2;
    
    % If simple or cockriging with one constraint, means are added back
    if itype<3
        t = t' + ones(m,1)*avg;
    else
        t = t';
    end
    x0s = t;
    
    % Calculation of cokriging variances
    
    s = kron(ones(m,1),sv);
    t = zeros(p,m);
    t(:) = diag(l'*k0);
    s = s-t';
end


        
        
    