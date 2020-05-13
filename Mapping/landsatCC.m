function [CC, NDVI] = landsatCC (Lat,Lon,buffer)
[A,Refmat,Box]=geotiffread('NWDIM.tif');
%proj= geotiffinfo('NWDIM.tif');% Projections of the tiff images of what we have.
mapshow(A,Refmat);% visualize the georectified image

%Lat=36.072;
%Lon=-97.212;

[X, Y] = latlon2utm (Lat, Lon, '14N');% changing the Coordiantes into UTM projection
%%to see whether the lat long that we converts in meters correctly falls on the image or not.
%mapshow (X,Y,'DisplayType','point','MarkerSize',20);

%imcrop% when this function is operated, it crop small subset image from large image as well as
rect = [5176.5 3347.5 435 344];% give xmin,ymin,width,height
Asub=imcrop(A,rect);% provides the pixel value of subset image.
%% d.Read Band 7, Band 5 and Band 4 from the subset image
Band5 = im2single(Asub(:,:,2));
Band4 = im2single(Asub(:,:,3));

%% Calculation of NWDI and NDVI
NDVI = (Band5 - Band4)./(Band5+Band4);% map the vegetation content in agricultural lands.

