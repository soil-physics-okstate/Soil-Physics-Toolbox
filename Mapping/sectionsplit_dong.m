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
% Modified 23-Jun-2015 by Geano Dong

%% Load sections map
sections_filename = ['Paynesections_UTM_14',filesep,'sections.shp'];
sections = shaperead(sections_filename);

%% Split sections. 

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

%% Iterate over each polygon

% Pre-allocate variables
centroid_qs_X = []; centroid_qs_Y = []; % X and Y coordinates of quarter section centroid
qsec_X = []; qsec_Y = []; % X and Y coordinates of the boundaries of each quarter section
C = []; % ramdom color for testing, can be ignored

for i=1:size(poly_X,1)
    % allocate coordinates for shorter name
    X = poly_X{i}; 
    Y = poly_Y{i}; 
    
    % exclude triangles 
    % (using '=' because of the 1st and the last elements are always the same)
    if length(X)<=4
        qsec_X = cat(1,qsec_X,X);      
        qsec_Y = cat(1,qsec_Y,Y);
        C = cat(1,C,rand(1,length(X)));
        centroid_qs_X = cat(1,centroid_qs_X,mean(X));
        centroid_qs_Y = cat(1,centroid_qs_Y,mean(Y));
        continue
    end
    
    % generate polygon edge vectors
    vecX = X(2:end)-X(1:end-1);
    vecY = Y(2:end)-Y(1:end-1);
    
    % calculate inner products for every two adjacent vectors in one polygon
    inprd = -vecX(1:end-1).*vecX(2:end)-vecY(1:end-1).*vecY(2:end);
    inprd = [-vecX(end).*vecX(1)-vecY(end).*vecY(1),inprd];
    
    % calculate the lengths of edges
    edges = sqrt(vecX(1:end).^2+vecY(1:end).^2);
    
    % calculate the cosine of each interior angle in one polygon
    cosIN = inprd./edges./[edges(end),edges(1:end-1)];
    
    % determine if a vertex is a right angle, approximately
    vert = cosIN<cosd(80) & cosIN>cosd(100);
    
    % extract the vertices with right angles
    vecXv = X(1:end-1);
    vecXv = vecXv(vert);
    vecYv = Y(1:end-1);
    vecYv = vecYv(vert);
    
%     % exclude triangles
%     if numel(vecXv)<4
%         qsec_X = cat(1,qsec_X,X);      
%         qsec_Y = cat(1,qsec_Y,Y);
%         C = cat(1,C,rand(1,length(X)));
%         centroid_qs_X = cat(1,centroid_qs_X,mean(X));
%         centroid_qs_Y = cat(1,centroid_qs_Y,mean(Y));
%         continue
%     end
    
    % calculate the centroid of this section
    sectcentroidX = mean(vecXv);
    sectcentroidY = mean(vecYv);
    
    % determine the centroid coordinates for each 'long edge'
    % generate vectors between adjacent vertices
    vecXsub = vecXv(2:end) - vecXv(1:end-1);
    vecXsub = [vecXsub,vecXv(1)-vecXv(end)];
    vecYsub = vecYv(2:end) - vecYv(1:end-1);
    vecYsub = [vecYsub,vecYv(1)-vecYv(end)];
    
    newedges = sqrt(vecXsub(1:end).^2+vecYsub(1:end).^2);
    if numel(newedges(newedges<820))>=2
        qsec_X = cat(1,qsec_X,X);      
        qsec_Y = cat(1,qsec_Y,Y);
        centroid_qs_X = cat(1,centroid_qs_X,mean(X));
        centroid_qs_Y = cat(1,centroid_qs_Y,mean(Y));
        C = cat(1,C,rand(1,length(X)));
        continue
    end
    
    % calculate the half-lengths of each 'long edge'
    halflength = 0.5*sqrt(vecXsub(1:end).^2+vecYsub(1:end).^2);
    
    % preallocation
%     Xnew = zeros(1,length(X)+length(vecXv));
%     Ynew = zeros(1,length(Y)+length(vecYv));
    Xnew = zeros(1,2*length(vecXv));
    Ynew = zeros(1,2*length(vecYv));
    n = 1;
    switches = [];
    order = [];
    for j = 1:length(X) % loop of polygon vertices
%         Xnew(n) = X(j);
%         Ynew(n) = Y(j);
%         n = n+1;
        for k = 1:length(vecXv) % loop of section corners
            if X(j) == vecXv(k)
                Xnew(n) = X(j);
                Ynew(n) = Y(j);
                n = n+1;
                if k < length(vecXv)
                    Xnew(n) = 0.5*(vecXv(k)+vecXv(k+1));
                    Ynew(n) = 0.5*(vecYv(k)+vecYv(k+1));
                    order = [order,n];
                    n = n+1;
                else
                    Xnew(n) = 0.5*(vecXv(k)+vecXv(1));
                    Ynew(n) = 0.5*(vecYv(k)+vecYv(1));
                    order = [order,n];
                    n = n+1;
                end
%                 if edges(k) > 800
%                      Xnew(n) = vecXv(k) + vecX(k)./edges(k).*halflength(k);
%                      Ynew(n) = vecYv(k) + vecY(k)./edges(k).*halflength(k);
%                     order = [order,n];
%                     n = n+1;
%                 elseif k < length(vecXv)
%                     Xnew(n) = vecXv(k+1)-vecX(k)./edges(k).*halflength(k);
%                     Ynew(n) = vecYv(k+1)-vecY(k)./edges(k).*halflength(k);
%                     switches = [switches, n];
%                     order = [order,n];
%                     n = n+1;
%                 else
%                     Xnew(n) = vecXv(1)-vecX(k)./edges(k).*halflength(k);
%                     Ynew(n) = vecYv(1)-vecY(k)./edges(k).*halflength(k);
%                     switches = [switches,n];
%                     order = [order,n];
%                     n = n+1;
%                 end
            end    
        end
    end
    % switch point positions
%     for m = switches
%         if m == length(Xnew)
%             Xnew = [Xnew(1),Xnew(m),Xnew(2:end-1)];
%             Ynew = [Ynew(1),Ynew(m),Ynew(2:end-1)];
%             order(order==m)=2;
%         else
%             tempX = Xnew(m);
%             Xnew(m) = Xnew(m+1);
%             Xnew(m+1) = tempX;
%             tempY = Ynew(m);
%             Ynew(m) = Ynew(m+1);
%             Ynew(m+1) = tempY;
%             order(order==m)=m+1;
%         end
%         clear tempX tempY
%     end
   
   % subloop construction
   qsectX = cell(length(order),1);
   qsectY = cell(length(order),1);
   qsectcentX = zeros(length(order),1);
   qsectcentY = zeros(length(order),1);
   c =  cell(length(order),1);
   for jj = 1:length(order)
       if jj == length(order) && order(1)~=1 
           qsectX{jj} = [sectcentroidX,Xnew(order(jj):end),Xnew(1:order(1))];
           qsectY{jj} = [sectcentroidY,Ynew(order(jj):end),Ynew(1:order(1))];
           qsectcentX(jj) = mean(qsectX{jj});
           qsectcentY(jj) = mean(qsectY{jj});
           c{jj} = rand(1,length(qsectX{jj}));
       else
           qsectX{jj} = [sectcentroidX,Xnew(order(jj):order(jj+1))];
           qsectY{jj} = [sectcentroidY,Ynew(order(jj):order(jj+1))];
           qsectcentX(jj) = mean(qsectX{jj});
           qsectcentY(jj) = mean(qsectY{jj});
           c{jj} = rand(1,length(qsectX{jj}));
       end
   end
   
   qsec_X = cat(1,qsec_X,qsectX);      
   qsec_Y = cat(1,qsec_Y,qsectY);
   centroid_qs_X = cat(1,centroid_qs_X,qsectcentX);
   centroid_qs_Y = cat(1,centroid_qs_Y,qsectcentY);
   C = cat(1,C,c);
end

%% Creating and saving maps

% Save KML file
% Overlaying KML files on Google maps require working with the World Geodesic
% System 1984 (WSG84) datum instead of the North American Datum 1983
% (NAD83). Working with NAD83 the quarter sections are shifted about 100
% meter to the north not matching well Google maps. So, conversion from UTM
% to WGS84 is done below using a function called utm2wgs that I obtained
% from the File Exchange (sorry I di not code this one).

% Convert X and Y into lat and lon (using WGS84 not NAD83)
qsec_lat = [];
qsec_lon = [];
for i=1:length(qsec_X)
    [lat,lon]=utm2deg(qsec_X{i}',qsec_Y{i}',repmat('14 S',length(qsec_Y{i}'),1));
    qsec_lat = cat(1,qsec_lat,{lat});
    qsec_lon = cat(1,qsec_lon,{lon});
end


% Set data in correct format (cell to structure arrays) and plot to ensure
% everything looks fine.
geometry = repmat({'Line'},length(qsec_X),1);
qs_geo = [qsec_lat qsec_lon];
S_geo = cell2struct([qs_geo geometry],{'Lat','Lon','Geometry'},2);

% disp('Saving KML file')
% kmlwrite('test.kml',S_geo) % It takes a while, although the file for%% Plot vector map
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
% Payne county is about 1MB (not sure why).

%% test fill.m
figure
hold on
for i = 1:size(qsec_X,1)
    fill(qsec_X{i},qsec_Y{i},C{i});
end
hold off
axis equal
