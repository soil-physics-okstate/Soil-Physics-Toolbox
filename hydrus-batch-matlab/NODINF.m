% NODINF.m
% Created by Derek Groenendyk
% 5/8/2012
% This a class that accesses the NOD_INF.OUT file

classdef NODINF < handle
    properties
        expFileLocation
        fid
        times
        numTimes
        firstNodes
        step
        averageHead
        lines
        numLayers
    end
      
    methods
        function nodinf = NODINF(directory)
            nodinf.expFileLocation = directory;
            nodinf.update()
            nodinf.readData()            
        end
        
        function update(nodinf)
            nodinf.fid = fopen([nodinf.expFileLocation,'\NOD_INF.OUT'],'r+');
            C = textscan(nodinf.fid, '%s', 'delimiter', '','whitespace', '');
            nodinf.lines = C{1};
            fclose(nodinf.fid);
        end

        function writeData(nodinf)
            frewind(nodinf.fid)
            for i=1:length(nodinf.lines)
                wline = nodinf.lines;
                fprintf(nodinf.fid,'%s\n',wline);
            end
            fclose(nodinf.fid);
        end  
        
        function readData(nodinf)
            nodinf.update()
            for i=1:length(nodinf.lines)
                if length(nodinf.lines{i}) > 1
                    sLine = textscan(nodinf.lines{i}, '%s');
                    if strcmp('Time:',sLine{1}{1})
                        startLine = i ; 
                    elseif strcmp('end',sLine{1}{1})
                        endLine = i;
                        break                       
                    end                    
                end           
            end
            
            nodinf.numLayers = endLine - startLine - 3;
            nodinf.step = endLine - startLine + 1;
            nodinf.numTimes = (length(nodinf.lines)-5)/nodinf.step;
            
            for i=1:nodinf.numTimes
               nodinf.times(end+1) = 6+(i-1)*nodinf.step;               
            end
            nodinf.firstNodes = nodinf.times + 3;

        end
              
        function allData = getAllData(nodinf)
            allData = zeros(nodinf.numTimes,nodinf.numLayers,5);        
            for ii=1:nodinf.numTimes
                for jj=1:nodinf.numLayers
                    index = nodinf.firstNodes(ii)+jj-1;
                    sLine = textscan(nodinf.lines{index}, '%s');
                    for kk=1:5
                        allData(ii,jj,kk) = str2double(sLine{1}{kk+2});
                    end
                end
            end
        end
        
    end % methods
end % classdef
