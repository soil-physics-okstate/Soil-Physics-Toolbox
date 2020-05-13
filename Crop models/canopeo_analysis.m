function [image_classified,CC] = canopeo_analysis(image_uint8,RG_ratio,BG_ratio,noise)
%%processing_core is the core function that classifies green canopy cover in Canopeo.
%
% Inputs
%       image_uint8 = 8 bits image.
%       RG_ratio = Typically a value between 0.9 and 1.1.
%       BG_ratio = Typically a value between 0.9 and 1.1.
%       noise = Try values 1, 10, 100, 1000. The higher the more connected
%               pixels (areas) are removed such as small weeds).
%
% Outputs
%       BW2 = binary image where true values (white pixels) represent green
%             canopy cover and false values (black pixels) represent background.
%       canopy = % Ground canopy cover.
%
% Supported image extension files: .jpg, .JPEG, .tiff, .TIF 
%
% Version 1.0 beta version
% 
% Andres Patrignani, 13 Feb 2014 9:28 AM.

green_tolerance = 20;
image_single = single(image_uint8); % Convert file from uint8 to single in order to use less memory.
RGratios = image_single(:,:,1)./image_single(:,:,2); % Ratio of red to green.
BGratios = image_single(:,:,3)./ image_single(:,:,2); % Ratio of blue to green.
tempClassified = (RGratios < RG_ratio) & (BGratios < BG_ratio) & (image_single(:,:,2) > green_tolerance); % Parameters for green color
image_classified = bwareaopen(tempClassified, noise); % Remove small areas.
CC = round((100*sum(sum(image_classified))/numel(image_classified))*100)/100; % Calculate percent canopy cover of the image.
