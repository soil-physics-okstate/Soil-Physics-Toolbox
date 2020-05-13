% OBSNODE.m
% Created by Derek Groenendyk
% 5/8/2012
% This a class that accesses the NOD_INF.OUT file

classdef OBSNODE < handle
    properties
        expFileLocation
        fid
        times
        numTimes
        numObs
        firstNodes
        step
        averageHead
        lines
        numLayers
    end
      
    methods
        function obsnode = OBSNODE(directory)
            obsnode.expFileLocation = directory;
            obsnode.update()
%             obsnode.readData()            
        end
        
        function update(obsnode)
            obsnode.fid = fopen([obsnode.expFileLocation,'\Obs_Node.out'],'r+');
            C = textscan(obsnode.fid, '%s', 'delimiter', '','whitespace', '');
            obsnode.lines = C{1};
            fclose(obsnode.fid);
        end

        function writeData(obsnode)
            frewind(obsnode.fid)
            for i=1:length(obsnode.lines)
                wline = obsnode.lines;
                fprintf(obsnode.fid,'%s\n',wline);
            end
            fclose(obsnode.fid);
        end  
        
        function readData(obsnode)
            obsnode.update()
            for i=1:length(obsnode.lines)
                if length(obsnode.lines{i}) > 1
                    sLine = textscan(obsnode.lines{i}, '%s');
                    if strcmp('Time:',sLine{1}{1})
                        startLine = i ; 
                    elseif strcmp('end',sLine{1}{1})
                        endLine = i;
                        break                       
                    end                    
                end           
            end
            
            obsnode.numLayers = endLine - startLine - 3;
            obsnode.step = endLine - startLine + 1;
            obsnode.numTimes = (length(obsnode.lines)-5)/obsnode.step;
            
            for i=1:obsnode.numTimes
               obsnode.times(end+1) = 6+(i-1)*obsnode.step;               
            end
            obsnode.firstNodes = obsnode.times + 3;

        end
              
        function obsData = getObsData(obsnode)
            % count num
            sLine = textscan(obsnode.lines{6}, '%s');
            numObs = length(sLine{1})/2;
            numTimes = length(obsnode.lines) - 8;
            
%             obsNode.numObs
%             obsNode.numTimes
            
            obsData = zeros(numTimes,numObs,2);
            for ii=1:numTimes
                sLine = textscan(obsnode.lines{ii+7}, '%s');
                for jj=1:numObs
                    obsData(ii,jj,1) = str2double(sLine{1}{2+(jj-1)*3});
                    obsData(ii,jj,2) = str2double(sLine{1}{3+(jj-1)*3});
                end
            end
        end
        
    end % methods
end % classdef
