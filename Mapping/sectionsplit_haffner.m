%%This script finds the centroids for quarter section of lands in Oklahoma.
%
% Matt, please feel free to modify this script. I will probably be working
% on it little more.
%
% Author Andres Patrignani - 17-Jun-2015

%% Load sections map
sections_filename = ['Paynesections',filesep,'sections.shp'];
sections = shaperead('OKsections_geo.shp');
% figure,mapshow(sections); hold on

%% RIVERS

% Load map of rivers
rivers_filename = ['Paynerivers_UTM_14',filesep,'rivers.shp'];
rivers = shaperead(rivers_filename);

% Project coordinates of rivers
X_rivers = [];
Y_rivers = [];
bufwidth = 0.005; %meters
for i=1:length(rivers)
    [lat_cell,lon_cell] = utm2latlon(rivers(i).X,rivers(i).Y,'14S');
    [latb,lonb] = bufferm(lat_cell,lon_cell,bufwidth);
    [X_cell,Y_cell] = latlon2utm(latb,lonb,'14S');
    rivers_area(i,1) = polyarea(X_cell(~isnan(X_cell)),Y_cell(~isnan(Y_cell))) / 1e6;
    X_rivers = cat(1,X_rivers,{X_cell});
    Y_rivers = cat(1,Y_rivers,{Y_cell});
    mapshow(X_cell,Y_cell);hold on;
end
rivers_area_threshold = 0.01; % 0.1 km2
rivers_area_idx = rivers_area > rivers_area_threshold;
X_rivers = X_rivers(rivers_area_idx);
Y_rivers = Y_rivers(rivers_area_idx);

%% LAKES

% Load map of lakes
lakes_filename = ['Paynelakes_UTM_14',filesep,'owrb_lakes.shp'];
lakes = shaperead(lakes_filename);

% Project coordinates of rivers
X_lakes = [];
Y_lakes = [];
for i=1:length(lakes)
    X_cell = lakes(i).X;
    Y_cell = lakes(i).Y;
    lakes_area(i,1) = polyarea(X_cell(~isnan(X_cell)),Y_cell(~isnan(Y_cell))) / 1e6;
    X_lakes = cat(1,X_lakes,{X_cell});
    Y_lakes = cat(1,Y_lakes,{Y_cell});
end

lakes_area_threshold = 0.1; % 0.1 km2
lakes_area_idx = lakes_area > lakes_area_threshold;
X_lakes = X_lakes(lakes_area_idx);
Y_lakes = Y_lakes(lakes_area_idx);

%% Split sections
X_poly = [];
Y_poly = [];
for i=1:size(sections,1)
    [X_cells,Y_cells] = polysplit(sections(i).X,sections(i).Y);
    X_poly = cat(1,X_poly,X_cells);
    Y_poly = cat(1,Y_poly,Y_cells);
end

%% Iterate over each section

% Pre-allocate variables
centroids_X = [];
centroids_Y = [];

for i=1:size(X_poly,1)

    % Find minimum of X and Y coordinates.
    X = X_poly{i}; % allocate X coordinate for shorter name
    min_X = min(X);
    max_X = max(X);
    span_X = (max_X-min_X); % East-West section length

    Y = Y_poly{i}; % allocate Y coordinate for shorter name
    min_Y = min(Y);
    max_Y = max(Y);
    span_Y = (max_Y-min_Y); % North-South section length

    % Find area of polygon
    section_area = span_X*span_Y / 1e6; %Calculate area in km squared
    
    % Find midpoints (currently not used)
    mid_X = span_X/2;
    mid_Y = span_Y/2;
    
    % Find 4, 2, or 1 centroid/s based on the area of the section.
    if section_area > 2 % A value of 2 km2 was arbitrarily chosen because 
                        % one section has a total of 2.56 km2 

        % Find centroids for large sections
        edges_X = linspace(min_X,max_X,5);
        centroids_X = cat(1, centroids_X, repmat(edges_X([2,4]),[1,2])');

        edges_Y = linspace(min_Y,max_Y,5);
        centroids_Y = cat(1, centroids_Y,[edges_Y([2,4]), edges_Y([4,2])]');
        
    elseif section_area > 1 && span_X > span_Y %smaller than 2 km2 and predominantly horizontal section
        
        % Find centroids
        edges_X = linspace(min_X,max_X,5); % 5 edges of the segmentes because is the longest side
        centroids_X = cat(1, centroids_X, edges_X([2,4])');

        edges_Y = linspace(min_Y,max_Y,3); % 3 edges of the segmentes because is the shortest side
        centroids_Y = cat(1, centroids_Y, [edges_Y(2), edges_Y(2)]');
        
     elseif section_area > 1 && span_X < span_Y
         
        % Find centroids
        edges_X = linspace(min_X,max_X,3); % 5 edges of the segmentes because is the longest side
        centroids_X = cat(1, centroids_X, [edges_X(2) edges_X(2)]');

        edges_Y = linspace(min_Y,max_Y,5); % 3 edges of the segmentes because is the shortest side
        centroids_Y = cat(1, centroids_Y, edges_Y([2,4])');
    
    elseif section_area > 0.2 && section_area < 1
        
        % Find centroids
        %edges_X = linspace(min_X,max_X,3);
        centroids_X = cat(1, centroids_X, mean(X));

        %edges_Y = linspace(min_Y,max_Y,3);
        centroids_Y = cat(1, centroids_Y, mean(Y));
    else
        continue
    end
end

%% Check if centroids are within water bodies

%RIVERS
for i=1:sum(rivers_area_idx)
    IN_rivers(i,:) = inpolygon(centroids_X,centroids_Y,X_rivers{i},Y_rivers{i});
end
IN_sum_rivers = sum(IN_rivers);

% LAKES
for i=1:sum(lakes_area_idx)
    IN_lakes(i,:) = inpolygon(centroids_X,centroids_Y,X_lakes{i},Y_lakes{i});
end
IN_sum_lakes = sum(IN_lakes);

IN_idx = IN_sum_rivers > 0 | IN_sum_lakes > 0 ;
centroids_X = centroids_X(~IN_idx);
centroids_Y = centroids_Y(~IN_idx);

%% Save (uncomment to run)
% centroids = [centroids_X(:),centroids_Y(:)];
% filename = 'centroids.csv';
% csvwrite(filename,centroids);

%% Plot (uncomment to run)
% Beware that these maps have lots points and may take a while to render.
% SAVE documents before plotting any map!
% opengl hardware
scatter(centroids_X(:),centroids_Y(:),'.r'); hold on
% mapshow(rivers); hold on;
mapshow(lakes)
