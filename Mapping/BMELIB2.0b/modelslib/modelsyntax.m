function modelsyntax

% modelsyntax               - Syntaxical help for using and creating models (Jan 01,2001)
%
% modelsyntax provides a help file to explain the syntax of the covariance
% models in BMELIB works.
%
% SYNTAX :
%
% help modelsyntax;
%
% or, more simply, you may just type:
%
% modelsyntax;
%
%
% CONTENT OF HELP FILE:
%
% 1- Basic Covariance models
% --------------------------
% 
% A covariance model in BMELIB allows to compute the covariance matrix 
% based on the distances between two sets of coordinates. The covariance
% model is speficied by its name covmodel, and its parameters covparam.
%
% Example 1.1:
%  Consider an exponential covariance for a 1-Dimentional domain with
%  sill cc=1 and range aa=1. A plot of this covariance is created as follow
%   cc=1;
%   aa=1;
%   covmodel='exponentialC';
%   covparam=[cc aa];
%   c1=[0];
%   c2=[0:0.1:2]';
%   [K]=coord2K(c1,c2,covmodel,covparam);
%   figure;
%   plot(c2,K);
%   xlabel('Spatial lag r');
%   ylabel('Covariance C(r)');
%
% Note that in BMELIB, the covariance is calculated using the function
%   [K]=coord2K(c1,c2,covmodel,covparam);
% This function calculates the distance matrix between the set of 
% coordinates c1 and c2 as follow
%   [D]=coord2dist(c1,c2);
% and then it computes the covariance matrix by using covmodel as follow
%   [K]=eval([covname '(D,covparam)']);
%
% Example 1.2:
%  Consider the same exponential covariance as above. The covariance
%  matrix can be directly computed from the distance matrix as follow
%   [D]=coord2dist(c1,c2);
%   [K]=exponentialC(D,covparam);
%
% It is worthwhile emphasizing as seen in these examples that 
% the covariance model supported by BMELIB are calculated based on 
% the distance between two sets of points. It is conceivable that the
% a user may change the covariance model by directly providing her/his
% own coord2K function. Such a function should then directly calculate
% the covariance values based on the set of coordinates c1 and c2.
% However as long as working with distances is acceptable, one should
% take advantage of the BMELIB built in covariance models. The 
% covariance models included with BMELIB are the following:
%
% nuggetC                   - nugget effect covariance model
% sphericalC                - spherical covariance model
% exponentialC              - exponential covariance model
% gaussianC                 - Gaussian covariance model
% holecosC                  - cosinusoidal hole effect covariance model
% holesinC                  - sinusoidal hole effect covariance model
%
% For each of the above covariance model the syntax to calculate 
% a covariance matrix is covname(D,covparam), where covname is the 
% name of the covariance function, D is the matrix of distances, and
% covparam is a vector of parameter. This first element of the 
% parameter vector, i.e. covparam(1), is the sill of the covariance
% model. In the case of the spherical, exponential and gaussian
% covariance models, covparam (2) is the range of the covariance, i.e.
% the distance to reach 95% of the sill value.
%
%
% 2 - space/time covariance models
% --------------------------------
%
% A coordinate set c1 is expressed as a n by nd matrix, where 
% n is the number of point and nd is the dimension of the spatial
% domain. For example in the 2D domain, the origin is given by
%   c1=[0 0];
% This concept is extended to the space time domain by specifying
% an additional column which holds the time. For e.g. the origin
% of a 2D spatial domain at time t=1 is given by
%   c1=[0 0 1];
% Whether a point is spatial only or space/time is determined 
% by the covariance model. In BMELIB the space/time covariance
% models comes in two flavors: separable, and non-separable.
% A space/time separable covariance model is coded as 
% 'covmodelS/covmodelT', where covmodelS refers to the spatial 
% covariance function and covmodelT refers to the temporal 
% covariance function. A non-separable covariance model is coded 
% as 'covmodelST'
%
% Example 2.1
%  Consider an space/time separable covariance model with a sill
%  of cc=1, a spatial gaussian component with spatial range as=1
%  in a 2D spatial domain, and a temporal exponential domain with 
%  a temporal range of at=2. A space/time color plot of this 
%  covariance is created as follow
%   cc=1;
%   as=1;
%   at=2;
%   covmodel='gaussianC/exponentialC';
%   covparam=[cc as at];
%   c1=[0 0 0];
%   sg=(0:0.1:1.5)';
%   nsg=length(sg);
%   tg=(0:0.1:3)';
%   ntg=length(tg);
%   c2=[kron(ones(ntg,1),sg) zeros(nsg*ntg,1) kron(tg,ones(nsg,1))];
%   [K]=coord2K(c1,c2,covmodel,covparam);
%   sgMat=reshape(c2(:,1),nsg,ntg);
%   tgMat=reshape(c2(:,end),nsg,ntg);
%   Kmat=reshape(K,nsg,ntg);
%   figure;
%   pcolor(sgMat,tgMat,Kmat);
%   colorbar;
%   shading interp;
%   title('Covariance C(r,t)');
%   xlabel('Spatial lag r');
%   ylabel('Temporal lag t');
%
% In this example the '/' character in the name of the covariance
% indicates that the model is space/time separable. As a result
% the last row of coordinates is taken as time, while the other rows
% are for the spatial position. Hence the two sets of space/time 
% coordinates c1 and c2, it is possible to calculate a matrix of 
% spatial distances Ds, as well as a matrix of temporal distances Dt.
% In other words, the calculation of the covariance matrix as follow
%   [K]=coord2K(c1,c2,'gaussianC/exponentialC',[cc aa at])
% is equivalent to calculating the spatial and time distance matrices,
% and then calling the covariance in space and time as follow
%   [Ds]=coord2dist(c1(:,1:end-1),c2(:,1:end-1));
%   [Dt]=coord2dist(c1(:,end),c2(:,end));
%   [K]=cc*gaussianC(Ds,[1 aa]).*exponentialC(Dt,[1 at]);
%
% Example 2.2
%  A space/time non-separable gaussian covariance model with a sill
%  of cc=1, and a space/time range of ast=1, where the space/time
%  distance is defined as follow
%    space/time distance = spatial distance + rst * temporal distance,
%  where rst=0.5 is a space/time metric, is given by the following 
%  covariance model:
%   cc=1;
%   ast=1;
%   rst=0.5;
%   covmodel='gaussianCST';
%   covparam=[cc ast rst];
%
% In this example the 'ST' characters in the name of the covariance
% indicates that the model is space/time separable. In this case 
% a spacial distance matrix Ds and a temporal distance matrix Dt
% are also calculated, but the covariance matrix is obtained as
%   [K]=gaussianCST(Ds,Dt,covparam);
%
%
% 3 - Nested covariance model
% --------------------------------
%
% When more than one covariance structure has to
% be specified, covmodel is a cell array where each
% cell is a string that contains the name of a model.
% Accordingly, covparam is a cell array where each cell
% contains the parameters of the corresponding covariance
% model.
%
% Example 3.1
%  To add a nugget effect with variance 0.2 to the 
%  space/time separable covariance model described in Example 2.1
%  we just use
%   cNugget=0.2;
%   cc=1;
%   as=1;
%   at=2;
%   covmodel={'nuggetC/nuggetC','gaussianC/exponentialC'};
%   covparam={ [cNugget], [cc as at] };
%
%
% 4 - Cross-Covariance model for several variables
% ------------------------------------------------
%
% It is possible to specify an additional index for the coordinates
% c1 and c2 taking integer values from 1 to nv. This index
% specifies which variable is known at each coordinate.
% In that case, c1 and c2 are cell arrays, where the first
% cell is the matrix of coordinates and the second cell is
% the vector of index values. Each covparam vector is then a
% cell array where the first cell is a symmetric covariance
% matrix for that model, and the second cell array is a
% vector that contains the other parameters of the models.
% The nhmax and nsmax variables are then vectors of nv
% elements, where sum(nsmax)<=20.
%
% Example 4.1
%  Consider two space/time variables X(s,t) and Y(s,t) in a
%  spatial domain of dimension 2, i.e. s=[s1 s2].
%  Consider the separable covariance model described in Example 3.1,
%  which has a nugget component nested with a exponentialC/gaussianC
%  separable component having a spatial range of as=1 and temporal 
%  range of at=2. Let's assume that the sills (cross-variances) of 
%  the nugget component for the variables X and Y is given by 
%    cNugget=[0.4 0.1; ...
%             0.1 0.1];
%  and the sill (cross-variances) for the exponentialC/gaussianC
%  component is given by
%    cc=[1  0.5; ... 
%        0.5 2];
%  A space/time color plot of the space/time cross-covariance 
%  between X and Y is created as follow
%   cNugget=[0.4 0.1;0.1 0.1];
%   cc=[1  0.5; 0.5 2];
%   as=1;
%   at=2;
%   covmodel={'nuggetC/nuggetC','gaussianC/exponentialC'};
%   covparam={ {cNugget,[]}, {cc, [as at]} };
%   c1=[0 0 0];
%   index1=[1];
%   p1={c1,index1};
%   sg=(0:0.1:1.5)';
%   nsg=length(sg);
%   tg=(0:0.1:3)';
%   ntg=length(tg);
%   c2=[kron(ones(ntg,1),sg) zeros(nsg*ntg,1) kron(tg,ones(nsg,1))];
%   n2=length(c2);
%   index2=2*ones(n2,1);
%   p2={c2,index2};
%   [K]=coord2K(p1,p2,covmodel,covparam);
%   sgMat=reshape(c2(:,1),nsg,ntg);
%   tgMat=reshape(c2(:,end),nsg,ntg);
%   Kmat=reshape(K,nsg,ntg);
%   figure;
%   pcolor(sgMat,tgMat,Kmat);
%   colorbar;
%   shading interp;
%   title('Cross covariance C_{X,Y}(r,t)');
%   xlabel('Spatial lag r');
%   ylabel('Temporal lag t');

help modelsyntax