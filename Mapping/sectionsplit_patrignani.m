%%This script:
%             1) Finds the centroids for each section
%             2) Divides each section into quarter sections
%             3) Finds centroids for each quarter section.
%             4) Saves KML, raster, and vector maps.
%
% The code in this document was tested using the shapefile with the section
% for Payne county, OK, but it is expected to work reasonably well at the
% US level.
%
% Feel free to modify or correct this script, just keep the green light on :)
%
% All maps should have datum NAD83 (ellipsoid GRS80) and projected UTM 14 coorindate system.
% For google maps we need to save maps using WGS84 datum (with ellipsoid
% WGS 84), but transformations are done before saving the maps (see bottom of
% the script for more details).

% RUN this script by setting your cwd in: \Soil Physics Toolbox\Mapping
%
% Author Andres Patrignani - 17-Jun-2015
% Modified 22-Jun-2015 (Added comments and polished code)

%% Load sections map
sections_filename = ['Paynesections_UTM_14',filesep,'sections.shp'];
%sections_filename = ['OKsections_UTM_14',filesep,'sections.shp'];
sections = shaperead(sections_filename);

% Split sections. 

% Some sections are composed of two or more polygons (e.g. sections near
% rivers or county edges). Separating the polygons allow to treat them
% separate. Polygons are divided at NaN values in the polygon sequence
% using the Matlab built-in function called polysplit.

poly_X = [];
poly_Y = [];
for i=1:size(sections,1)
    [X_cells,Y_cells] = polysplit(sections(i).X,sections(i).Y);

    % Collect splitted polygons
    poly_X = cat(1,poly_X,X_cells);
    poly_Y = cat(1,poly_Y,Y_cells);
end

poly_X = cellfun(@(x) [x,NaN],poly_X,'UniformOutput',false);
poly_Y = cellfun(@(x) [x,NaN],poly_Y,'UniformOutput',false);

% Iterate over each polygon

% Pre-allocate variables
centroid_qs_X = []; centroid_qs_Y = []; % X and Y coordinates of quarter section centroid
qsec_X = []; qsec_Y = []; % X and Y coordinates of the boundaries of each quarter section
dx=[]; dy=[]; % Approximate span in the X and Y direction for each polygon. 
              % Not very meaningful of non-squared polygons.
xy11 = []; % Variable to store the lower and left coordinate of each polygon.
section_area = nan(size(poly_X,1),1);

for i=1:size(poly_X,1)

    disp([num2str(i/size(poly_X,1)*100),' %']) % Time consuming task
    
    % Find convoluted polygon
    X = poly_X{i}; % allocate X coordinate for shorter name
    Y = poly_Y{i}; % allocate Y coordinate for shorter name
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
%     K = convhull(X,Y,'simplify',true);
%     X = X(K);
%     Y = Y(K);
    
    % Find minimum of X coordinates.
    min_X = min(X);
    max_X = max(X);
    dx = cat(1,dx,max_X-min_X); % East-West section length

    % Find minimum of Y coordinates.
    min_Y = min(Y);
    max_Y = max(Y);
    dy = cat(1,dy,max_Y-min_Y); % North-South section length
    xy11 = cat(1,xy11,[min_X min_Y]);
    
    [~,idxpoint] = unique([X;Y]','rows','stable');
    
    % Find area of polygon
    section_area(i,1) = dx(i)*dy(i) / 1e6; % Calculate area in km squared
       
    % An arbitrary value was chosen because one section has a total of 2.56 km2 
    if section_area(i,1)>2.5 && length(idxpoint)>=4  

        % Find centroid for full section
        edges_X = linspace(min_X,max_X,5);
        %centroid_s_X = edges_X(3);
        centroid_s_X = nanmean(X);
        edges_Y = linspace(min_Y,max_Y,5);
        %centroid_s_Y = edges_Y(3);
        centroid_s_Y = nanmean(Y);
        
        % Find centroids for quarter sections
        centroid_qs_X = cat(1, centroid_qs_X, repmat(edges_X([2,4]),[1,2])');        
        centroid_qs_Y = cat(1, centroid_qs_Y,[edges_Y([2,4]), edges_Y([4,2])]');
        
        % Initialize steps to split full section into quarter sections
        % Normalize polygon points relative to the X and Y mean
        % coordinates to devide the section in four cuadrants. This step
        % helps finding the four corners of the full section.
        Xn = X-centroid_s_X;
        Yn = Y-centroid_s_Y;
        
        % Calculate distance from full section centroid
        d = sqrt(Xn.^2+Yn.^2); 
        
        % Eliminate repeated points that are used to close the polygon
        [~,idxpoint] = unique([X;Y]','rows','stable');
        d = d(idxpoint);
        
        % Since some points that are used to create the polygon are close 
        % to the corners, and because some of the polygons are not
        % perfectly squared, some of these points could have greater
        % distances to the full section centroid than some of the other
        % corners. To fix this I decided to grab n points with the greatest
        % distance from the full section. n should be >= 4 (we typically need four
        % corners for sections with area > 2.5 km squared). This entire
        % script is to account for a small portion of cases, meaning that
        % mostly the 4 corners are being correctly selected with the above
        % procedure, but occasionally one of the corners appear in 5th
        % place. Setting n=5 allows manipulation and extraction of all 4
        % corners.
        % 4
        X2 = [X(idxpoint) X(idxpoint)];
        Y2 = [Y(idxpoint) Y(idxpoint)];
        N = length(d);
        angles = [];
        for j=2:N+1
            angles(j) = angle_rad( [X2(j-1);Y2(j-1)],[X2(j);Y2(j)],[X2(j+1);Y2(j+1)]);
        end
        angles(1) = angles(end); 
        angles(end) = [];
        [~,idxcorner] = sort(abs(angles),'descend');
        corner_Xn = Xn(idxcorner(1:4));
        corner_Yn = Yn(idxcorner(1:4));
        
        % Re-organize corners. Again, in most cases the points that
        % describe the polygon are correctly organized counter clockwise 
        % starting from the top left (NW) corner of the polygon, but in
        % some cases they are not, and so we need re-organize them.
        % The first step consists on finding in which quadrant each corner
        % is located.
        
        % Find cuadrant for each corner.
        idx_I = corner_Xn<0 & corner_Yn>0;
        idx_II = corner_Xn>0 & corner_Yn>0;
        idx_III = corner_Xn>0 & corner_Yn<0;
        idx_IV = corner_Xn<0 & corner_Yn<0;
        
        % Re-organize corners
        corner_X = [corner_Xn(idx_I) corner_Xn(idx_II) corner_Xn(idx_III) corner_Xn(idx_IV)]+centroid_s_X;
        corner_Y = [corner_Yn(idx_I) corner_Yn(idx_II) corner_Yn(idx_III) corner_Yn(idx_IV)]+centroid_s_Y;
        
        % Useful when debugging
        % mapshow(X,Y);hold on
        % scatter(corner_X(1),corner_Y(1),'sr','filled');hold on
        % scatter(corner_X(2:4),corner_Y(2:4),'sb','filled');hold on
        
        % With corners in the right place, now we need to create a
        % four-point polygon of the quarter section. Because the shape of
        % the full section is not a perfect square creating the quarter
        % section polygon is not trivial. Previously we re-organized
        % counter clockwise the corners (the corners of the full section), 
        % now we will need to repeat the same step with all four points of
        % the polygon describing the quarter section. NaNs are added at the
        % end to let Matlab know that the polygon ends up (allows for
        % storing all quarter section polygons together)
        
        qsec1_X = [corner_X(1) mean(corner_X([1,2])) centroid_s_X mean(corner_X([1,4])) corner_X(1)];
        qsec1_Y = [corner_Y(1) mean(corner_Y([1,2])) centroid_s_Y mean(corner_Y([1,4])) corner_Y(1)];
        [qsec1_X,qsec1_Y] = poly2cw(qsec1_X,qsec1_Y);
          
        qsec2_X = [corner_X(2) mean(corner_X([1,2])) centroid_s_X mean(corner_X([2,3])) corner_X(2)];
        qsec2_Y = [corner_Y(2) mean(corner_Y([1,2])) centroid_s_Y mean(corner_Y([2,3])) corner_Y(2)];
        [qsec2_X,qsec2_Y] = poly2cw(qsec2_X,qsec2_Y);

        qsec3_X = [corner_X(3) mean(corner_X([2,3])) centroid_s_X mean(corner_X([3,4])) corner_X(3)];
        qsec3_Y = [corner_Y(3) mean(corner_Y([2,3])) centroid_s_Y mean(corner_Y([3,4])) corner_Y(3)];
        [qsec3_X,qsec3_Y] = poly2cw(qsec3_X,qsec3_Y);

        qsec4_X = [corner_X(4) mean(corner_X([3,4])) centroid_s_X mean(corner_X([1,4])) corner_X(4)];
        qsec4_Y = [corner_Y(4) mean(corner_Y([3,4])) centroid_s_Y mean(corner_Y([1,4])) corner_Y(4)];
        [qsec4_X,qsec4_Y] = poly2cw(qsec4_X,qsec4_Y);
        
        % Combine all points describing the quarter section polygons.
        qsec_X = cat(1,qsec_X,{qsec1_X;qsec2_X;qsec3_X;qsec4_X});
        qsec_Y = cat(1,qsec_Y,{qsec1_Y;qsec2_Y;qsec3_Y;qsec4_Y});
        
    else    
        % No division into quarter sections for full section with area
        % lower than the threshold.
        % Find and store quarter section centroids
        centroid_qs_X = cat(1, centroid_qs_X, nanmean(X));
        centroid_qs_Y = cat(1, centroid_qs_Y, nanmean(Y));

        [qsec_X_temp,qsec_Y_temp] = poly2cw(X,Y);
        [qsec_X_temp,qsec_Y_temp] = closePolygonParts(qsec_X_temp,qsec_Y_temp);
        
        % Store full section polygon boundaries as in the shapefile.
        qsec_X = cat(1,qsec_X,qsec_X_temp);
        qsec_Y = cat(1,qsec_Y,qsec_Y_temp);

    end
end


%% Convert X and Y into lat and lon (using WGS84 not NAD83)
qsec_lat = [];
qsec_lon = [];
qcentroid_lat = [];
qcentroid_lon = [];
for i=1:length(qsec_X)
    [lat,lon] = utm2deg(qsec_X{i}',qsec_Y{i}',repmat('14 S',length(qsec_Y{i}'),1));
    qsec_lat = cat(1,qsec_lat,{lat});
    qsec_lon = cat(1,qsec_lon,{lon});
    
    [lat2,lon2] = utm2deg(centroid_qs_X(i),centroid_qs_Y(i),'14 S');
    qcentroid_lat = cat(1,qcentroid_lat,{lat2});
    qcentroid_lon = cat(1,qcentroid_lon,{lon2});    
end

%% Save in JSON format

% Set data in correct format (cell to structure arrays) and plot to ensure
% everything looks fine.
geometry = repmat({'Line'},length(qsec_X),1);
qs_geo = [qsec_lat qsec_lon qcentroid_lat qcentroid_lon];
S_geo = cell2struct([qs_geo geometry],{'Qsec_lat','Qsec_lon','Qcentroid_lat','Qcentroid_lon','Geometry'},2);

jsonobj = savejson('payne_qsection',S_geo,'payne_qsection');

% Save KML file
% Overlaying KML files on Google maps require working with the World Geodesic
% System 1984 (WSG84) datum instead of the North American Datum 1983
% (NAD83). Working with NAD83 the quarter sections are shifted about 100
% meter to the north not matching well Google maps. So, conversion from UTM
% to WGS84 is done below using a function called utm2wgs that I obtained
% from the File Exchange.

% disp('Saving KML file')
% kmlwrite('test.kml',S_geo) % It takes a while, although the file for
% Payne county is about 1MB (not sure why).

%% Plot vector map
% If plotting state level maps, beware that may take a while for your 
% computer to render the maps. SAVE your documents.
% opengl software % Can help rendering the maps.

qs_map = [qsec_X qsec_Y];
S_map = cell2struct([qs_map geometry],{'X','Y','Geometry'},2);
figure,
mapshow(sections); hold on
mapshow(S_map); hold on
scatter(centroid_qs_X,centroid_qs_Y,'.r'); hold on
axis equal

%% Plot raster map
% Raster maps can have some distortion when working with non-uniform grids
% Unfortunately, the reference matrix R only accepts one pixel spacing 
% between consecutive centroids.

x11 = min(centroid_qs_X);
y11 = min(centroid_qs_Y);
dx = median(dx(section_area>2.5)/2); % Median pixel distance in the X direction.
dy = median(dy(section_area>2.5)/2); % Median pixel distance in the Y direction.
R = makerefmat(x11,y11,dx,dy);
[row,col] = map2pix(R,centroid_qs_X,centroid_qs_Y);

% Output from map2pix is not an integer in case you were expecting the same
% I was. Matlab returns a float in case operations within the pixel are
% desired (it is easier to round than to recover the exact position within
% the pixel).
row = round(row);
col = round(col);

% Crete raster image
A = zeros(max(row),max(col),3);
% Sorry, A(row(i),col(i),2) = 1; was not working properly for me. Hence, a
% lowely loop.
for i=1:length(row)
    A(row(i),col(i),2) = 1;
end
% A = uint8(A);
figure
mapshow(A,R); hold on;
mapshow(S_map); hold on
% geotiffname = 'test.tif';
% geotiffwrite(geotiffname,A,R);

%% Save CSV with full section centroids
% centroids = [centroids_X(:),centroids_Y(:)];
% filename = 'centroids.csv';
% csvwrite(filename,centroids);
