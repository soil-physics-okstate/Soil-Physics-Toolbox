function plotZsBME(estvalues,mapwhat,mapfun,transfun,celloptions,vecoptions,harddata,softdata,mask,areas)

% plotZsBME                 - Create maps showing the spatial distribution of a BME estimate
%
% creates a map showing the spatial distribution of the Space/Time Random Field (S/TRF)
% Z(s,t) at a given time t=tk.  The map shows the distribution of the estimated 
% values Z(s,t) using either contour lines or a color shading, as well as possibly the
% hard and soft data points.
%
%   In general the estimated values available and the hard and soft data may correspond 
% to a transform X(s,t) of the actual variable Z(s,t).  For example in the case of 
% a log-transform, then X(s,t) is the log-transform of Z(s,t), i.e. X(s,t)=log(Z(s,t).
% In this case the estimated values and the hard and soft data for X(s,t) are provided
% as input, and the function can plot maps of X(s,t), or maps of Z(s,t) if the 
% back-transform function is specified (e.g. the back transform of the log-transform
% is Z(s,t)=exp(X(s,t)) )
%
%   The map generated can be any of the following, depending on the mapwhat input variable
% 1) Map of the estimated value (either X(s,t) or it's back-transform Z(s,t))
% 2) Map of the length of the Confidence Interval (CI) for either X(s,t) or Z(s,t)
%    at a specified confidence level
% 3) Map of any quantile of X(s,t) or Z(s,t) at any confidence probability
% 4) Map of the normalized estimation error variance of X(s,t), (i.e. the estimation error
%    variance of X(s,t) divided by the variance of X(s,t))
% 5) Maps delineating non attainment area at a given confidence level for Z(s,t), i.e.
%    areas not attaining Z(s,t)<Zlimit within some confidence level 
%
% SYNTAX :
%
% plotZsBME(estvalues,mapwhat,mapfun,transfun,celloptions,vecoptions,harddata,softdata,mask,areas)
%
% INPUT :
%
% estvalues={sk,XkBME,XkBMEv,tk,varianceX,varname,varunit,xlab,ylab}; or [] or filename
%
%                           estvalues contains the BME mean estimate of the (possibly transformed)
%                           Spatial/Temporal Random Field (S/TRF) X(s,t).  It is either a cell array   
%                           with the input values shown above, or an empty set, or a character string 
%                           with the name of a file.  When estvalues is empty the estimated value is  
%                           missing and it is not plotted.  Otherwise the input values are provided by 
%                           the cell array or from the file.  If an input value is not provided by the 
%                           cell array or the file then it's default value is used. The input values
%                           of estvalues are the following:
%   sk          nk x 2      spatial coordinates of the estimation point
%   XkBME       nk x 1      BME estimated value of the S/TRF X(s,t) at the estimation points
%   XkBMEv      nk x 1      estimation error variance of the S/TRF X(s,t) at the estimation points
%   tk          scalar      optional input value with the time of the monitoring event, used to write the
%                           title. if tk=[] then this input value is not used. The default value is []
%   varianceX   scalar      optional input value with variance of the data X. if varianceX=NaN then the 
%                           variance is calculated from the XkBME data. The default value is NaN
%   varname     char str    optional input value with the name of the variable. If varname='' then the default
%                           value is used.  The default value is 'Z'
%   varunit     char str    optional input value with the unit of Z. If varunit='' then the unit is not used.
%                           The default value is ''
%   xlab        char str    character string with the label of the x-axis.  The default is 's1'
%   ylab        char str    character string with the label of the y-axis.  The default is 's2'
% 
%
% mapwhat       scalar      indicates what to map, as follow
%                           1 map of the estimated value (either X(s,t) or it's back-transform Z(s,t))
%                           2 map of the length of the Confidence Interval (CI) for either X(s,t) or Z(s,t)
%                             at a specified confidence level
%                           3 map of any quantile of X(s,t) or Z(s,t) at any confidence probability
%                           4 map of the normalized estimation error variance of X(s,t), the  nestimation error
%                             variance of X(s,t) divided by the variance of X(s,t)
%                           5 maps delineating non attainment area at a given confidence level for Z(s,t), i.e.
%                             areas not attaining Z(s,t)<Zlimit within some confidence level 
%                           The default value is 1
%
%
% mapfun        scalar      indicates which plotting function to use to create the maps
%                           1 uses contour.m to create contours
%                           2 uses pcolor to create color checherboard maps
%                           3 uses contourf.m to create filled contours
%                           4 uses mesh to create a 3D mesh surface
%                           5 uses surf to create a 3D surface
%                           11 uses contour.m to create contours, without labels
%                           13 uses contourf.m to create filled contours, without labels
%                           The default value is 2
%
%
% transfun={invtransfun,transname,plottrans}; or []
%
%                           transfun specifies the transformation used to obtain X(s,t) from the
%                           original variable Z(s,t).  transfun is either a cell array of input values
%                           as shown above, or an empty set.  If transfun is empty then this indicate
%                           that there was no transformation used.  If transfun is defined, then a
%                           transformation was used, and the estimated values and hard and soft data
%                           provided are that of X(s,t)=g(Z(s,t)) where g(.) is the transformation
%                           function.  In this case transfun has the following input values :
%   invtransfun  char str   character string with the name of the inverse of the transformation
%                           function.  For example is X(s,t) is the log-transform of Z(s,t), then
%                           invtransfun='exp';  
%   transname    char name  character string specifying the name of the transformation function. This
%                           name is only used for the title of the figure.   For example if X(s,t) is 
%                           the log-transform of Z(s,t), then transname='log';
%   plottrans    scalar     indicates what scale to use to plot plot X(s,t) 
%                           (plottrans=1) plot the back-transform Z(s,t)=invtransfun(X(s,t)) .
%                           (plottrans=2) plot X(s,t) 
%                           (plottrans=3) plot Z(s,t)=invtransfun(X(s,t)) on a log scale
%
%
% celloptions={ax,cax,cmap,contval,griddatamethod}; or []
%
%                           celloptions specifies some optional parameters controling the making of the
%                           maps. celloptions is either a cell array of input values
%                           as shown above, or an empty set.  If celloptions is empty then default values 
%                           are used. Otherwise celloptions defines the following optional parameters :
%   ax        4x1 or empty  ax=[xmin xmax ymin ymax] specifies the spatial area to view.  If ax is empty
%                           then the viewing area is set to comprise all the estimation and data points.
%                           The default value is []
%   cax       2x1 or empty  cax=[minZ maxZ] (or [minX maxX]) specifies the lower and upper bound cuttof
%                           values used in the colorbar of the color maps of the variable mapped.
%                           If cax is empty then the min and max of the estimated values (and of the data
%                           if defined) are used.
%                           The default value is []
%   cmap      n x 3         A lookup table defining the colormap of the color maps.  See help colormap
%                           for additional information.
%   contval   1 x m         A vector of values used to draw the contour lines if mapfun=1. See 
%                           help contour for additional information.
%                           By default the contour.m automatically generates contour values
%   griddatamethod          char string defining the option used in the griddata.m function.
%                           See help griddata for the possible values of this option.
%                           The default is griddatamethod='linear';
%
% vecoptions=[conflevel,nonattainval,nxpix,nypix,showestpoints,showlegend]; or []
%
%                           vecoptions specifies some optional parameters controling the making of the
%                           maps. vecoptions is either a cell array of input values
%                           as shown above, or an empty set.  If vecoptions is empty then default values 
%                           are used. Otherwise vecoptions defines the following optional parameters :
%   conflevel      scalar   value indicating a confidence level as follow
%                           if mapwhat=2 or 4,  confidence level of confidence intervals
%                           if mapwhat=3, confidence level of the quantile to plot
%                           if mapwhat=5, confidence level of the attainment areas
%                           The default is 0.68 if mapwhat=2 or 4, 0.8 if mapwhat=3 and 0.5 if mapwhat=3.
%   nonattainval   scalar   limit value of used to delineate non-attainment area when mapwhat=5. 
%                           The non-attainment areas are areas not attaining Z<nonattainval withing the 
%                           specified confidence level (i.e. where P(Z<nonattainval)>conflevel is false)
%                           Note that the value specified is for Z(s,t).  When using transformed values
%                           X(s,t) then estvalues, harddata and softdata are values for X(s,t), but 
%                           nonattainval is the non attainment value for the backtransform Z(s,t).
%                           The default value is the average of available estimation values or hard and
%                           soft data.
%   nxpix          scalar   integer specifying the number of pixel in the x direction used to create the
%                           colored maps. The default value is 200
%   nypix          scalar   integer specifying the number of pixel in the y direction used to create the
%                           colored maps. The default value is 160
%   showestpoints  scalar   optional integer indicating whether to show the estimation points:
%                           1 to show estimation points, 0 ow
%                           The default is 0.
%   showestpoints  scalar   optional integer indicating whether to show a legend of the points:
%                           1 to show the legend, 0 ow
%                           The default is 1.
%
%
% harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim}; or [] or sh or filename
%
%                           harddata contains the hard data of the (possibly transformed)
%                           Spatial/Temporal Random Field (S/TRF) X(s,t).  It is either a cell array   
%                           with the input values shown above, or an empty set, or a nhx2 numerical
%                           variable, or a character string with the name of a file.
%                           When harddata is empty the hard data is missing and it is not plotted.
%                           When harddata is a nhx2 numerical matrix is defines sh, the spatial coordinates
%                           of the hard data points, and only the location of the points are plotted.
%                           Otherwise the input values are provided by the cell array or from the file,
%                           allowing to plot the hard data points with markerplot or colorplot.
%                           If an input value is not provided by the cell array or the file then it's default
%                           value is used. The input values of estvalues are the following:
%   sh            nh x 2    spatial coordinates of the hard data points
%   zh            nh x 1    values of S/TRF X(s,t) at the hard data points.  If zh is empty then only
%                           the location of the points can be shown.  The default is zh=[]
%   shlegend      char str  character string to define hard data in the legend. 
%                           Default is shlegend='hard data points'
%   shplotfun     scalar    integer specifying the plotting function used to display the hard data, as follow
%                           1) to plot just the location of the hard data points
%                           2) to use markerplot (see help markerplot)
%                           3) to use colorplot (see help colorplot)
%                           4) to use poleplot (see help poleplot)
%                           The default is 1
%   shProperties  cell      cell array defining 'Property' in markerplot.m or colorplot.m
%   shValues      cell      cell array defining 'Value' in markerplot.m or colorplot.m
%   shzrange      2 x 1     cell array defining 'zrange' in markerplot.m or colorplot.m
%   shsizelim     2 x 1     cell array defining 'sizelim' in markerplot.m or colorplot.m
%
% 
% softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim}; or [] or ss or filename
%
%                           Same as harddata, but for the soft data.
%
%
% mask={maskcontour,maskfillcolor,masklinetype}; or empty or filename
%
%                           mask contains a maskcontour around which the area is filled
%                           with a maskfillcolor. mask is either a cell array   
%                           with the input values shown above, or an empty set, or a nx2 numerical
%                           variable, or a character string with the name of a file.
%                           When mask is empty the maskcontour is missing and it is not plotted.
%                           When mask is a nx2 numerical matrix is defines the maskcontour and
%                           a default color is used. Otherwise the input values are provided by the
%                           cell array or from the file. If an input value is not provided by the 
%                           cell array or the file then it's default value is used. The input 
%                           values of estvalues are the following:
%
%   maskcontour    n x 2    A nx2 matrix of points defining the outer contour of the mask
%   maskfillcolor  char     a character defining a color, see the plot function.
%                           The default value is 'w' (white)
%   masklinetype   char str a character string defining the linetype of mask contour
%                           The default value is 'k' (black line)
%
% areas={areascontour,areaslinetype,areaslegend}; or [] or filename
%
%
%                           areas delineates the boundaries of areas. areas is either a cell array   
%                           with the input values shown above, or an empty set, or a character 
%                           string with the name of a file.
%                           When areas is empty the area boundaries are missing and are not plotted.
%                           Otherwise the input values are provided by the cell array or from 
%                           the file. If an input value is not provided by the cell array 
%                           or the file then it's default value is used. The input 
%                           values of areas are the following:
%
%   areascontour  cell      a cell array of areas boundaries.  Each cell element is a nx2 matrix 
%                           of points defining the outer contour of that area
%   masklinetype  char str  a character string defining the linetype used to draw the area contours
%                           The default value is ':k' (dotted black line) for mapfun=1
%                           and 'k' (black line) for mapfun=2
%   areaslegend   char str  a character string witht he name of the area boundaries to use for the
%                           legend. The default value is areaslegend='area boundaries'
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default input parameters
%
if nargin<1, estvalues='estvalues.mat', end;
if nargin<2, mapwhat=1; end;
if nargin<3, mapfun=2; end;
if nargin<4, transfun=[]; end;
if nargin<5, celloptions=[]; end;
if nargin<6, vecoptions=[]; end;
if nargin<7, harddata=[]; end;
if nargin<8, softdata=[]; end;
if nargin<9, mask=[]; end;
if nargin<10, areas=[]; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input parameters for errors and extract values from each input parameters
%

%
%  Get the estimation values from estvalues={sk,XkBME,XkBMEv,tk,varname,varunit,varianceX,xlab,ylab}
%
if isempty(estvalues)
  sk=[];
  XkBME=[];
  XkBMEv=[];
  tk=[];
  nk=0;
  varianceX=NaN;
  varname='';
  varunit='';
  xlab='';
  ylab='';
elseif ischar(estvalues)
  load(estvalues);
  if ~exist('sk','var') | ~exist('XkBME','var') | ~exist('XkBMEv','var') 
    error(sprintf('The file %s must contain the variables sk, XkBME, XkBMEv',estvalues));
  end
  if ~exist('tk','var'), tk=[]; end
  if ~exist('varianceX','var'), varianceX=NaN; end
  if ~exist('varname','var'), varname=''; end
  if ~exist('varunit','var'), varunit=''; end
  if ~exist('xlab','var'), xlab='s1'; end
  if ~exist('ylab','var'), ylab='s2'; end
  nk=size(sk,1);
  dim=size(sk,2);
elseif iscell(estvalues), 
  if ~isnumeric(estvalues{1}), error('estvalues{1} must be numeric'); end;
  sk=estvalues{1};
  nk=size(sk,1);
  dim=size(sk,2);
  if dim~=2, error('sk=estvalues{1} must have two columns'); end;
  if length(estvalues)<2
    error('estvalues must be a cell array of at least two cells'); 
  else
    if ~isnumeric(estvalues{2}), error('estvalues{2} must be numeric'); end;
    XkBME=estvalues{2};
    if size(XkBME,1)~=nk, error('XkBME=estvalues{2} must have the same number of lines as sk=estvalues{1}'); end;
    if length(XkBME)>0 & size(XkBME,2)~=1, error('XkBME=estvalues{2} must have only one column'); end;
  end
  if length(estvalues)<3
    XkBMEv=[];
  else
    if ~isnumeric(estvalues{3}), error('estvalues{3} must be numeric'); end;
    XkBMEv=estvalues{3};
    if size(XkBMEv,1)~=nk, error('XkBMEv=estvalues{3} must have the same number of lines as sk=estvalues{1}'); end;
    if length(XkBMEv)>0 & size(XkBMEv,2)~=1, error('XkBMEv=estvalues{3} must have only one column'); end;
  end
  if length(estvalues)<4
    tk=[];
  else 
    if ~isnumeric(estvalues{4}), error('estvalues{4} must be numeric'); end;
    tk=estvalues{4};
    if length(tk)>1, error('tk=estvalues{4} must be a scalar'); end;
  end
  if length(estvalues)<5,
    varianceX=NaN;
  else
    varianceX=estvalues{5};
  end
  if length(estvalues)<6,
    varname='Z';
  else
    if ~ischar(estvalues{6}), error('estvalues{6} must be a char string'); end;
    varname=estvalues{6};
  end
  if length(estvalues)<7,
    varunit='';
  else
    if ~ischar(estvalues{7}), error('estvalues{7} must be a char string'); end;
    varunit=estvalues{7};
  end
  if length(estvalues)<8,
    xlab='';
  else
    if ~ischar(estvalues{8}), error('estvalues{8} must be a char string'); end;
    xlab=estvalues{8};
  end
  if length(estvalues)<9,
    ylab='';
  else
    if ~ischar(estvalues{9}), error('estvalues{9} must be a char string'); end;
    ylab=estvalues{9};
  end
else
  error('estvalues must be either [], or a char string, or a cell array'); 
end
  

% get the integer mapfun the mapping function to use
switch mapfun
case {11,13}, 
  mapfun=mapfun-10;
  writeclabel=0;
otherwise, 
  writeclabel=1;
end


%
%  Get the transformation function from transfun={invtransfun,transname,plottrans};
%
if isempty(transfun)
  invtransfun='';
  transname='';
  plottrans=0;
elseif iscell(transfun),
  if ~ischar(transfun{1}), error('transfun{1} must be char string'); end;
  invtransfun=transfun{1};
  if length(transfun)<3, error('transfun must be a cell array of at least three cells'); end
  if ~ischar(transfun{2}), error('transfun{2} must be a char string'); end;
  transname=transfun{2};
  if ~isnumeric(transfun{3}) | length(transfun{3})~=1, error('transfun{3} must be a scalar'); end;
  plottrans=transfun{3};
  if plottrans~=1 & plottrans~=2 & plottrans~=3, error('transfun{3} must be equal to 1 or 2 or 3'); end
else
  error('transfun must be either [] or a cell array'); 
end
if isempty(invtransfun), 
  invtransfun='identity';
  transname='identity';
end


%
%  Get the options from celloptions={ax,cax,cmap,contval,griddatamethod};
%
if isempty(celloptions)
  ax=[];
  cax=[];
  cmap=[];
  contval=[];
  griddatamethod='linear';
elseif iscell(celloptions),
  if ~isnumeric(celloptions{1}), 
    error('celloptions{1} must be numeric');
  else
    ax=celloptions{1};
    if length(ax)>0,
      [n1 n2]=size(ax);
      if n1~=1 | n2~=4, error('celloptions{1} must be a 1 by 4 vector'); end;
    end
  end    
  if length(celloptions)<2, 
    cax=[];
  elseif ~isnumeric(celloptions{2}), 
    error('celloptions{2} must be numeric'); 
  else
    cax=celloptions{2};
    if length(cax)>0,
      [n1 n2]=size(cax);
      if n1~=1 | n2~=2, error('celloptions{2} must be a 1 by 2 vector'); end;
    end
  end;
  if length(celloptions)<3, 
    cmap=[];
  elseif ~isnumeric(celloptions{3}), 
    error('celloptions{3} must be numeric'); 
  else
    cmap=celloptions{3};
    [n1 n2]=size(cmap);
    if n1==0 | n2~=3, error('celloptions{3} must be a n by 3 matrix'); end;
  end;
  if length(celloptions)<4, 
    contval=[];
  elseif ~isnumeric(celloptions{4}), 
    error('celloptions{4} must be numeric'); 
  else
    contval=celloptions{4};
    if length(contval)>0,
      [n1 n2]=size(contval);
      if min([n1 n2])~=1, error('celloptions{4} must be a vector'); end;
    end
  end;
  if length(celloptions)<5, 
    griddatamethod='linear';
  elseif ~ischar(celloptions{5}), 
    error('celloptions{5} must be a char string'); 
  else
    griddatamethod=celloptions{5};
  end;
else
  error('celloptions must be either [] or a cell array'); 
end

%
%  Get the options from vecoptions=[conflevel,nonattainval,nxpix,nypix,showestpoints,showlegend];
%
if isempty(vecoptions)
  conflevel=NaN;
  nonattainval=NaN;
  nxpix=200;
  nypix=160;
  showestpoints=0;
  showlegend=1; 
elseif isnumeric(vecoptions),
  if length(vecoptions)<1 
    conflevel=NaN; 
  else
    conflevel=vecoptions(1);     
  end
  if length(vecoptions)<2, 
    nonattainval=NaN; 
  else
    nonattainval=vecoptions(2);     
  end
  if length(vecoptions)<3, 
    nxpix=100; 
  else
    nxpix=vecoptions(3);     
  end
  if length(vecoptions)<4, 
    nypix=80; 
  else
    nypix=vecoptions(4);     
  end
  if length(vecoptions)<5, 
    showestpoints=0; 
  else
    showestpoints=vecoptions(5);     
  end
  if length(vecoptions)<6, 
    showlegend=1; 
  else
    showlegend=vecoptions(6);     
  end
else
  error('vecoptions must be either [] or a numeric vector'); 
end

%
%  Get the hard data from
%      harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim};
%
if isempty(harddata)
  sh=[];
  nh=0;
  zh=[];
  shlegend='hard data points';
  shProperties={};
  shValues={};
  shplotfun=1;
  shzrange=[];
  shsizelim=[];
elseif isnumeric(harddata)
  sh=harddata;
  nh=size(sh,1);
  dim=size(sh,2);
  if dim~=2, error('sh=harddata must have two columns'); end;
  zh=[];
  shlegend='hard data points';
  shProperties={};
  shValues={};
  shplotfun=1;
  shzrange=[];
  shsizelim=[];
elseif ischar(harddata)
  load(harddata);
  if ~exist('sh','var') error(sprintf('The file %s must contain the variables sh',harddata)); end
  if ~exist('zh','var') zh=[]; end
  if ~exist('shlegend','var') shlegend='hard data points'; end
  if ~exist('shplotfun','var') shplotfun=1; end
  if ~exist('shProperties','var') shProperties=[]; end
  if ~exist('shValues','var') shValues=[]; end
  if ~exist('shzrange','var') shzrange=[]; end
  if ~exist('shsizelim','var') shsizelim=[]; end
  nh=size(sh,1);
elseif iscell(harddata), 
  if ~isnumeric(harddata{1}), error('harddata{1} must be numeric'); end;
  sh=harddata{1};
  nh=size(sh,1);
  dim=size(sh,2);
  if dim~=2, error('sh=harddata{1} must have two columns'); end;
  if length(harddata)<2, 
    zh=[];
  else
    if ~isnumeric(harddata{2}), error('harddata{2} must be numeric'); end;
    zh=harddata{2};
    if size(zh,1)~=nh, error('zh=harddata{2} must have the same number of lines as sh=harddata{1}'); end;
    if length(zh)>0 & size(zh,2)~=1, error('zh=harddata{2} must have only one column'); end;
  end
  if length(harddata)<3, 
    shlegend='hard data points';
  else
    if ~ischar(harddata{3}), error('harddata{3} must be a char string'); end;
    shlegend=harddata{3};
  end
  if length(harddata)<4,
    shplotfun=1;
  else
    if ~isnumeric(harddata{4}) | length(harddata{4})~=1, error('harddata{4} must be a scalar'); end;
    shplotfun=harddata{4};
  end
  if length(harddata)<5,
    shProperties={};
  else
    if ~iscell(harddata{5}), error('harddata{5} must be a cell array'); end;
    shProperties=harddata{5};
  end
  if length(harddata)<6,
    shValues={};
  else
    if ~iscell(harddata{6}), error('harddata{6} must be a cell array'); end;
    shValues=harddata{6};
  end
  if length(harddata)<7, 
    shzrange=[];
  elseif ~isnumeric(harddata{7}), 
    error('harddata{7} must be numeric'); 
  else
    shzrange=harddata{7};
    if length(shzrange)>0,
      [n1 n2]=size(shzrange);
      if n1~=1 | n2~=2, error('harddata{7} must be a 1 by 2 vector'); end;
    end
  end;
  if length(harddata)<8, 
    shsizelim=[];
  elseif ~isnumeric(harddata{8}), 
    error('harddata{8} must be numeric'); 
  else
    shsizelim=harddata{8};
    if length(shsizelim)>0,
      [n1 n2]=size(shsizelim);
      if n1~=1 | n2~=2, error('harddata{8} must be a 1 by 2 vector'); end;
    end
  end;
else
  error('harddata must be either [], or a nh by 2 matrix, or a char string, or a cell array'); 
end
  
%
%  Get the soft data from
%      softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim};
%
if isempty(softdata)
  ss=[];
  ns=0;
  zs=[];
  sslegend='soft data points';
  ssProperties={};
  ssValues={};
  ssplotfun=1;
  sszrange=[];
  sssizelim=[];
elseif isnumeric(softdata)
  ss=softdata;
  ns=size(ss,1);
  dim=size(ss,2);
  if dim~=2, error('ss=softdata must have two columns'); end;
  zs=[];
  sslegend='soft data points';
  ssProperties={};
  ssValues={};
  ssplotfun=1;
  sszrange=[];
  sssizelim=[];
elseif ischar(softdata)
  load(softdata);
  if ~exist('ss','var') error(sprintf('The file %s must contain the variables ss',softdata)); end
  if ~exist('zs','var') zs=[]; end
  if ~exist('sslegend','var') sslegend='soft data points'; end
  if ~exist('ssplotfun','var') ssplotfun=1; end
  if ~exist('ssProperties','var') ssProperties=[]; end
  if ~exist('ssValues','var') ssValues=[]; end
  if ~exist('sszrange','var') sszrange=[]; end
  if ~exist('sssizelim','var') sssizelim=[]; end
  ns=size(ss,1);
elseif iscell(softdata),
  if ~isnumeric(softdata{1}), error('softdata{1} must be numeric'); end;
  ss=softdata{1};
  ns=size(ss,1);
  dim=size(ss,2);
  if dim~=2, error('ss=softdata{1} must have two columns'); end;
  if length(softdata)<2, 
    zs=[];
  else
    if ~isnumeric(softdata{2}), error('softdata{2} must be numeric'); end;
    zs=softdata{2};
    if size(zs,1)~=ns, error('zs=softdata{2} must have the same number of lines as ss=softdata{1}'); end;
    if length(zs)>0 & size(zs,2)~=1, error('zs=softdata{2} must have only one column'); end;
  end
  if length(softdata)<3, 
    sslegend='soft data points';
  else
    if ~ischar(softdata{3}), error('softdata{3} must be a char string'); end;
    sslegend=softdata{3};
  end
  if length(softdata)<4,
    ssplotfun=1;
  else
    if ~isnumeric(softdata{4}) | length(softdata{4})~=1, error('softdata{6} must be a scalar'); end;
    ssplotfun=softdata{4};
  end
  if length(softdata)<5,
    ssProperties={};
  else
    if ~iscell(softdata{5}), error('softdata{5} must be a cell array'); end;
    ssProperties=softdata{5};
  end
  if length(softdata)<6,
    ssValues={};
  else
    if ~iscell(softdata{6}), error('softdata{6} must be a cell array'); end;
    ssValues=softdata{6};
  end
  if length(softdata)<7, 
    sszrange=[];
  elseif ~isnumeric(softdata{7}), 
    error('softdata{7} must be numeric'); 
  else
    sszrange=softdata{7};
    if length(sszrange)>0,
      [n1 n2]=size(sszrange);
      if n1~=1 | n2~=2, error('softdata{7} must be a 1 by 2 vector'); end;
    end
  end;
  if length(softdata)<8, 
    sssizelim=[];
  elseif ~isnumeric(softdata{8}), 
    error('softdata{8} must be numeric'); 
  else
    sssizelim=softdata{8};
    if length(shsizelim)>0,
      [n1 n2]=size(shsizelim);
      if n1~=1 | n2~=2, error('softdata{8} must be a 1 by 2 vector'); end;
    end
  end;
else
  error('softdata must be either [], or a ns by 2 matrix, or a char string, or a cell array'); 
end
  
if nk==0 & nh==0 & ns==0
  return;
end

%
%  Get the mask from mask={maskcontour,maskfillcolor,masklinetype};
%
if  isempty(mask)
  maskcontour=[];
  maskfillcolor=[];
elseif isnumeric(mask)
  maskcontour=mask;
  [n1 n2]=size(maskcontour);
  if n1<=3, error('mask must have at least three lines'); end
  if n2~=2, error('mask must have two columns'); end
  maskfillcolor='w';
  masklinetype='k';
elseif ischar(mask)
  load(mask);
  if ~exist('maskcontour'), mask=[]; end
  if ~exist('maskfillcolor'), maskfillcolor='w'; end
  if ~exist('masklinetype'), masklinetype='k'; end
elseif iscell(mask),
  if ~isnumeric(mask{1}), error('mask{1} must be numeric'); end;
  maskcontour=mask{1};
  [n1 n2]=size(maskcontour);
  if n1<=3, error('mask{1} must have at least three lines'); end
  if n2~=2, error('mask{1} must have two columns'); end
  if length(mask)<2, 
    maskfillcolor='w';
  elseif ~char(mask{2}), 
    error('mask{2} must be a char'); 
  else
    maskfillcolor=mask{2};
  end;
  if length(mask)<3, 
    masklinetype='k';
  elseif ~char(mask{3}), 
    error('mask{3} must be a char string'); 
  else
    masklinetype=mask{3};
  end;
else
  error('mask must be either [], or numeric, or a char string, or a cell array'); 
end

%
%  Get the mask from areas={areascontour,areaslinetype,areaslegend};
%
if  isempty(areas)
  areascontour=[];
  areaslinetype=[];
  areaslegend=[];
elseif ischar(areas)
  load(areas);
  if ~exist('areascontour') areascontour=[]; end
  if ~exist('areaslinetype') areaslinetype=':k'; end
  if ~exist('areaslegend') 
    if length(areascontour)<=1
      areaslegend='area boundary';
    else
      areaslegend='area boundaries';
    end
  end
else
  if ~iscell(areas), error('areas must be a cell array'); end
  areascontour=areas{1};
  if ~iscell(areas{1}), error('areas{1} must be a cell array'); end
  nareas=length(areascontour);
  for i=1:nareas
    if ~isnumeric(areascontour{i}), error('Each cell element of areas{1} must be numeric'); end
    [n1 n2]=size(areascontour{i});
    if n1<3, error('Each cell element of areas{1} must have at least three lines'); end
    if n2~=2, error('Each cell element of areas{1} must have two columns'); end
  end
  if length(areas)<2
    areaslinetype=':k';
  else
    areaslinetype=areas{2};
    if ~ischar(areas{2}), error('areas{2} must be a character string'); end
  end
  if length(areas)<3
    if length(areascontour)<=1
      areaslegend='area boundary';
    else
      areaslegend='area boundaries';
    end
  else
    areaslegend=areas{3};
    if ~ischar(areas{3}), error('areas{3} must be a character string'); end
  end
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign default values
%
clabelFontSize=8;
if isempty(xlab), xlab='s1'; end;
if isempty(ylab), ylab='s2'; end;
if isempty(griddatamethod), griddatamethod='linear'; end;       % linear interpolation of contours
if isempty(shsizelim), shsizelim=[5 20]; end;
if isempty(sssizelim), sssizelim=[5 20]; end;

if isempty(areaslinetype)
  if mapfun==2, 
    areaslinetype='k';
  else
    areaslinetype=':k';
  end
end

%
% Set new figure
%
figure;
hp=[];
hold on;

%
% Set default value for cmap
%
if isempty(cmap)
  cmap=hot;
  cmap=cmap(end:-1:1,:);
end

%
%    Set default value for conflevel
%
if isnan(conflevel),
  switch mapwhat
  case {2,4}, conflevel=0.68;     % confidence probability of confidence interval
  case 3, conflevel=0.80;         % confidence level of quantile
  case 5, conflevel=0.5;          % confidence level of non attainment zones
  end
end

%
%    Set default value for ax
%
if isempty(ax) 
  s=[];
  if nk>0, s=[s;sk]; end
  if nh>0, s=[s;sh]; end
  if ns>0, s=[s;ss]; end
  ax=[min(s(:,1)) max(s(:,1)) min(s(:,2)) max(s(:,2))];
end

%
%    Calculate the variable to plot Zk
%
switch(mapwhat)
case 1,
  if plottrans==0 | plottrans==2, Zk=XkBME;
  else, Zk=feval(invtransfun,XkBME);
  end;
  if plottrans==3, 
    idx=Zk>0;
    if sum(idx)==0, error('Need at least some positive values to plot on a log10-scale'); end;
    Zk(idx)=log10(Zk(idx)); 
    Zk(~idx)=floor(min(Zk(idx))-1);
  end
case 2
  confcoef=gaussinv((1+conflevel)/2,[0 1]);
  conflowbound=XkBME-confcoef*sqrt(XkBMEv);
  confupbound=XkBME+confcoef*sqrt(XkBMEv);
  if plottrans==0 | plottrans==2, 
    Zk=confupbound-conflowbound;
  else
    conflowbound=feval(invtransfun,conflowbound);
    confupbound=feval(invtransfun,confupbound);
    Zk=confupbound-conflowbound;
  end
  if plottrans==3, 
    idx=Zk>0;
    if sum(idx)==0, error('Need at least some positive values to plot on a log10-scale'); end;
    Zk(idx)=log10(Zk(idx)); 
    Zk(~idx)=floor(min(Zk(idx))-1);
  end
case 3
  confcoef=gaussinv(conflevel,[0 1]);
  quantileval=XkBME+confcoef*sqrt(XkBMEv);
  if plottrans==0 | plottrans==2, 
    Zk=quantileval;
  else
    Zk=feval(invtransfun,quantileval);
  end
  if plottrans==3, 
    idx=Zk>0;
    if sum(idx)==0, error('Need at least some positive values to plot on a log10-scale'); end;
    Zk(idx)=log10(Zk(idx)); 
    Zk(~idx)=floor(min(Zk(idx))-1);
  end
case 4
  if isnan(varianceX), varianceX=max(XkBMEv); end;
  Zk=XkBMEv/varianceX;
case 5
  confcoef=gaussinv(conflevel,[0 1]);
  quantileval=XkBME+confcoef*sqrt(XkBMEv);
  if plottrans==0, 
    Zk=quantileval;
  else
    Zk=feval(invtransfun,quantileval);
  end
otherwise
  error('mapwhat must be equal to 1, 2, 3, 4, or 5');
end

%
%    Set default value for nonattainval
%
if mapwhat==5 & isnan(nonattainval)
  nonattainval=mean(Zk(~isnan(Zk)));
end

%
%    Transform the zh and zs values if needed
%
if length(zh)>0 & (plottrans==1 | plottrans==3)
  zh=feval(invtransfun,zh);
  if plottrans==3, 
    idx=zh>0;
    if sum(idx)==0, error('Need at least some positive hard data to plot on a log10-scale'); end;
    zh(idx)=log10(zh(idx)); 
    zh(~idx)=floor(min(zh(idx))-1);
  end
end
if length(zs)>0 & (plottrans==1 | plottrans==3)
  zs=feval(invtransfun,zs);
  if plottrans==3, 
    idx=zs>0;
    if sum(idx)==0, error('Need at least some positive soft data to plot on a log10-scale'); end;
    zs(idx)=log10(zs(idx)); 
    zs(~idx)=floor(min(zs(idx))-1);
  end
end

%
%    Set default value for cax
%
if isempty(cax) 
  if nh>0 & ~isempty(shzrange)
    cax=shzrange;
  elseif ns>0 & ~isempty(sszrange) 
    cax=sszrange;
  elseif nk>0
    cax=[min(Zk) max(Zk)];
  elseif nh>0
    cax=[min(zh) max(zh)];
  elseif ns>0
    cax=[min(zs) max(zs)];
  end    
end
sszrange=cax;
shzrange=cax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Create the figure with the map  

if isempty(varname)
  varname='Z ';
else
  varname=[varname ' '];
end


if isempty(tk) | isnan(tk) 
  eventtime='';
else
  eventtime=sprintf('at time =%g ',tk);
end

disp(sprintf('Plotting map of %s %s',varname,eventtime));

%
% Construct the grid used for display and interpolate the values to plot
% on that grid
%
if nk>0
  dx1=diff(ax(1:2))/nxpix;
  dy1=diff(ax(3:4))/nypix;
  xg=[ax(1):dx1:ax(2)+dx1-eps];
  yg=[ax(3):dy1:ax(4)+dy1-eps];
  [xg yg]=meshgrid(xg,yg);
  Zg=griddata(sk(:,1),sk(:,2),Zk,xg,yg,griddatamethod);
  Zg=reshape(Zg,size(xg));
  maxZg=max(max(Zg));
end
  
%
% Plot the contour lines or pcolor map of the Zk data to plot
%
if mapwhat==5 & nk>0
  cax=[nonattainval-diff(cax)/10 nonattainval+diff(cax)/10];
  cmap=[1 1 1;1 0 0];
end
if ~isempty(cax)
  caxis(cax);
end
colormap(cmap); 
if nk>0
  switch mapfun
  case 1
    [hc hl]=contour(xg,yg,Zg,'-.k');
    if writeclabel==1 
      htclabel=clabel(hc,hl,'FontSize',clabelFontSize);
    end
  case 2
    pcolor(xg,yg,Zg);
    shading interp;
  case 3
    [hc hl]=contourf(xg,yg,Zg,'-.k');
    if writeclabel==1 
      htclabel=clabel(hc,hl,'FontSize',clabelFontSize);
    end
  case 4
    mesh(xg,yg,Zg);
  case 5
    surf(xg,yg,Zg);
  otherwise
    error('bad value for mapfun');
  end
end

%
% Change the scale on the contour labels if using the log10 scale
%
if exist('htclabel','var') & plottrans==3
  for i=1:length(htclabel)
    if isequal(get(htclabel(i),'Type'),'text'),
      yt=str2num(get(htclabel(i),'String'));
      yt=10.^yt;
      ytr=10.^floor(log10(yt))./100.*round(100*yt./(10.^floor(log10(yt))));
      str=sprintf('%3.2g',ytr);
      set(htclabel(i),'String',str);
    end
  end
end

%
% Show Hard Data  
%
if nh>0
  h=[];
  switch shplotfun
  case 1,  
    h=plot(sh(:,1),sh(:,2),'o');
    for j=1:length(shProperties),
      set(h,shProperties{j},shValues{j});
    end;
  case 2,  
    markerplot(sh,zh,shsizelim,shProperties,shValues,shzrange);
    if showlegend & ~isempty(shlegend), 
      h=plot(ax(1)-1,ax(3)-1,'o'); 
    end
  case 3,  
    colorplot(sh,zh,cmap,shProperties,shValues,shzrange); 
    if showlegend & ~isempty(shlegend), 
      h=plot(ax(1)-1,ax(3)-1,'o');
    end
  case 4,  
    poleplot(sh,zh,shProperties,shValues); 
  end
  if showlegend & ~isempty(shlegend) & ~isempty(h)
    for j=1:length(shProperties),
      set(h,shProperties{j},shValues{j});
    end;
    hp(length(hp)+1)=h;
    leg{length(hp)}=shlegend;
  end
end

%
% Show Soft Data
%
if ns>0
  h=[];
  switch ssplotfun
  case 1,  
    h=plot(ss(:,1),ss(:,2),'v');
    for j=1:length(ssProperties),
      set(h,ssProperties{j},ssValues{j});
    end;
  case 2,  
    markerplot(ss,zs,sssizelim,ssProperties,ssValues,sszrange);
    if showlegend & ~isempty(sslegend), 
      h=plot(ax(1)-1,ax(3)-1,'v'); 
    end
  case 3,  
    colorplot(ss,zs,cmap,ssProperties,ssValues,sszrange); 
    if showlegend & ~isempty(sslegend), 
      h=plot(ax(1)-1,ax(3)-1,'v');
    end
  case 4,  
    poleplot(ss,zs,ssProperty,ssValue); 
  end
  if showlegend & ~isempty(sslegend) & ~isempty(h)
    for j=1:length(ssProperties),
      set(h,ssProperties{j},ssValues{j});
    end;
    hp(length(hp)+1)=h;
    leg{length(hp)}=sslegend;
  end
end

%
% Show estimation points
%
if showestpoints & nk>0
  hp(length(hp)+1)=plot(sk(:,1),sk(:,2),'+k');
  leg{length(hp)}='Estimation points';
end;

%
% Write the title
%

if mapwhat~=4
  if plottrans==0 | plottrans==1 | isempty(transname)
    transname='';
  elseif plottrans==2
    transname=[transname '-transformed '];
  elseif plottrans==3
    transname=['log10-transformed '];
  end
else
  if plottrans==0 | isempty(transname)
    transname='';
  elseif plottrans==2
    transname=[transname '-transformed '];
  elseif plottrans==3
    transname=['log10-transformed '];
  end
end  

if ~isempty(varunit)
  varunit=['(' varunit ') '];
end

switch(mapwhat)
case 1,
  title(sprintf('Map of %s%s%s%s',transname,varname,varunit,eventtime));
case 2
  title(sprintf('Length of the %d%% CI of %s%s%s%s',conflevel*100,transname,varname,varunit,eventtime));
case 3
  title(sprintf('%d%% quantile of %s%s%s%s',conflevel*100,transname,varname,varunit,eventtime));
case 4
  title(sprintf('Normalized estimation error variance of %s%s%s%s',transname,varname,varunit,eventtime));
case 5
  title(sprintf('Areas not attaining an %d%% confidence level of %s<%g %s%s',...
    conflevel*100,varname,nonattainval,varunit,eventtime));
end

%
%  Fill in the color around the mask
%
if ~isempty(mask)
  axmask = [ min([maskcontour(:,1);ax(1)]) max([maskcontour(:,1);ax(2)])...
      min([maskcontour(:,2);ax(3)]) max([maskcontour(:,2);ax(4)]) ];
  [dummy,i]=min( abs( maskcontour(:,1) - axmask(1)));
  i=i(1);
  mask=[maskcontour(i:end,:);maskcontour(1:i-1,:)];
  idx=[find(sum(abs(diff(maskcontour,1,1)),2)~=0);size(maskcontour,1)];
  maskcontour=maskcontour(idx,:);
  mask=[maskcontour;maskcontour(1,:)];
  if maskcontour(2,2)>maskcontour(1,2), maskcontour=maskcontour(end:-1:1,:); end;
  fmask=[[axmask(1),maskcontour(1,2)];maskcontour;[axmask(1),maskcontour(1,2)];...
      [axmask(1) axmask(4)];[axmask(2) axmask(4)];[axmask(2) axmask(3)];...
      [axmask(1) axmask(3)];[axmask(1),maskcontour(1,2)]];
  h=fill(fmask(:,1),fmask(:,2),maskfillcolor);
  set(h,'EdgeColor',maskfillcolor);
  plot(maskcontour(:,1),maskcontour(:,2),masklinetype);
end

%
%  plot the areas
%
if ~isempty(areas)
  for i=1:length(areascontour)
    h=plot(areascontour{i}(:,1),areascontour{i}(:,2),areaslinetype);
  end
  hp(length(hp)+1)=h;
  leg{length(hp)}=areaslegend;
end

xlabel(xlab);
ylabel(ylab);
axis(ax);

%
% Put a colorbar and change its scale is using the log10 scale
%
if ( nk>0 & (mapfun==2 | mapfun==3) & mapwhat~=5) | (nh>0 & shplotfun==3) | (ns>0 & ssplotfun==3)
  hc=colorbar;
  if plottrans==3
    yt=get(hc,'YTick');
    yt=10.^yt;
    ytr=10.^floor(log10(yt))./100.*round(100*yt./(10.^floor(log10(yt))));
    for i=1:length(ytr)
      str=sprintf('%.5g',ytr(i));
      ytl(i,1:length(str))=str;
    end
    set(hc,'YTickLabel',ytl);  
  end
end
if ~isempty(hp) & showlegend
  legend(hp,leg);
end

function y=identity(x)
y=x;