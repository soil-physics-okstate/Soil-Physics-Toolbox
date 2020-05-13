%Created by TEO, November 2012
%First run "MakeDepthTable.m"
depthtable = dlmread('\tabular\DepthTableclay.txt',',',1,0); %this variable contains the map unit id and corresponding clay content for the surface horizon

%Set map limits in UTM coordinates
MapXLimit = [651108 668164];
MapYLimit = [3986604 3997264];
mapunits = shaperead('\spatial\soilmu_a_ok119.shp','BoundingBox', [MapXLimit' MapYLimit']);
%This loop adds the surface horizon clay content to each map unit polygon
for i = 1:length(depthtable)
    for j = 1:length(mapunits)
        if str2num(mapunits(j).MUKEY) == depthtable(i,1)
            mapunits(j).clay = depthtable(i,3);         
        end
    end
end

%Identify the maximum clay content for the area:
maxclay = max([mapunits.clay]);

%Create an autumn colormap for the area, and then use the flipud command to invert the matrix.
fall = flipud(autumn(numel(mapunits)));

%Make a symbol specification structure, a symbolspec, that assigns an
%autumn color to each polygon according to the clay content
clayColors = makesymbolspec('Polygon', {'clay', ...
   [0 maxclay], 'FaceColor', fall});

%Draw axes with the desired limits
axes('XLim',MapXLimit,'YLim',MapYLimit)

%Display the map.
mapshow(mapunits, 'DisplayType', 'polygon', ...
   'SymbolSpec', clayColors,'LineWidth',0.001)
set(gca, 'XTick', []);
set(gca, 'YTick', []);

%Create a colorbar.
caxis([0 maxclay])
colormap(fall)
colorbar
