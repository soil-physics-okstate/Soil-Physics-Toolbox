% Objective to retrieve Landsat-8 image of June,2013 using MATLAB tool.
% Calculate Normalized Difference Water Index(NDWI) and Normalized
% Diffrence Vegetation Index(NDVI). Landsat-8 is recently launced satellite
% images. It has altogether 11 bands.NDWI and NDVI are calculated from
% three bands i.e.Band 7,Band 5 and Band4.
% Imporvement made since version 1.0
% * Crop the image of Marena from large Landsat image.
% * Calculate NDWI and NDVI from the cropped image.
% * Shows the NDWI and NDVI maps of Marena
% Author = Sonisa Sharma
% Date = 8-20-2013
%% Steps that are followed
% a. Read the image with geotiffread command and then give
% the name A. Digital numbers from three bands are used as our input
% variables.
[b04,Refmat,Box]=geotiffread('b04.TIF');
[b05,Refmat1,Box1]=geotiffread('b05.TIF');
[b07,Refmat2,Box2]=geotiffread('b07.TIF');
% changing the uint 16 to uint 8 so that it can be converted to the value
% range from 0 to 255
b4 = uint8(b04./ 256);
b5 = uint8(b05./256);
b7 = uint8(b07./256);
% stack the image
compositeImage = cat(3,b7,b5,b4); 
imshow(compositeImage);
%
Lon=36.072;
Lat=-97.212;
% imcrop;
% rect = [654465 3994365 3660 4380];
% Asub = imcrop(A,rect);
% [X1,Y1] = map2pix(Refmat,rect(1),rect(2));% want to check
[X, Y] = latlon2utm (Lat, Lon, '14N');% changing the Coordiantes into UTM projection
%%to see whether the lat long that we converts in meters correctly falls on the image or not
mapshow(X,Y,'DisplayType','point','MarkerSize',30);

%% b.Provided Latitude and Longitude are converted into UTM zone 14 projection
Lon=36.072;
Lat=-97.212;
% imcrop;
% rect = [654465 3994365 3660 4380];
% Asub = imcrop(A,rect);
% [X1,Y1] = map2pix(Refmat,rect(1),rect(2));% want to check
[X, Y] = latlon2utm (Lat, Lon, '14N');% changing the Coordiantes into UTM projection
%%to see whether the lat long that we converts in meters correctly falls on the image or not.
mapshow(X,Y,'DisplayType','point','MarkerSize',20);
title('Georectified image')
%mapshow (X1,Y1,'DisplayType','point','MarkerSize',20);
%% c.Zoom in the image so that the area of interest i.e Marena can be drawn.
% image is used even though we see the place as in Marker using mapshow.
% The disadvantages of image is that it gave the image without projection
% but the advantages is it works well with imcrop to crop the areaof
% interest from large image.
image(A)% when this command is used, there is no marker but we know the place earlier using
% line 24.
title('Ungeorectified image');
imcrop% when this function is operated, it crop small subset image from large image as well as
rect = [5176.5 3347.5 435 344];% give xmin,ymin,width,height
Asub=imcrop(A,rect);% provides the pixel value of subset image.
%% d.Read Band 7, Band 5 and Band 4 from the subset image
% Band 7 is short wave infrared which has wavelength of 2.11-2.29
% micrometer, Band 5 is Near infrared  with wavelength of 0.85-0.88
% micrometer and Band 4 is the red band with wavelenth of 0.64-0.67
% micrometer. All the three bands are spatial resolution of 30m.
Band7 = im2single(Asub(:,:,1));% im2single converts the stack image into one band.
Band5 = im2single(Asub(:,:,2));
Band4 = im2single(Asub(:,:,3));
figure
subplot(2,2,1:2),imshow(Band7);% subplot shows all the images in one figure
% so that it is good for user to know what the image is he/she looking.
title('Short wave region')% provides the title of the image
subplot(2,2,3),imshow(Band5);% imshow provide the visualization of the image.
title('Near Infrared Band');
subplot(2,2,4),imshow(Band4);
title('Red Band');
%% Calculation of NWDI and NDVI
NWDI=(Band7 - Band5)./(Band7 + Band5);% map the water content in agricultural lands
NDVI = (Band5 - Band4)./(Band5+Band4);% map the vegetation content in agricultural lands.
subplot(2,1,1),imshow(NWDI);
title('Normalized Difference Water Index')
subplot(2,1,2),imshow(NDVI);
title('Normalized Difference Vegetation Index')
%% From Matlab central
img1 = imread('b04.TIF');
img2 = imread('b05.TIF');
img3 = imread('b07.TIF');

combined_img(:,1) = img1;
combined_img(:,2) = img2;
combined_img(:,3) = img3;
