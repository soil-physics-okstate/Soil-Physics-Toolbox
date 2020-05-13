% ATMOSPHIN.m
% Created by Derek Groenendyk
% created 1/27/2015
% This a class that accesses the ATMOSPH.IN file
% modified 1/27/2016 by Ben Paras

classdef ATMOSPHIN < handle
    properties
        expFileLocation
        fid
        C
    end
      
    methods        
        function atmosph = ATMOSPHIN(directory)
            atmosph.expFileLocation = directory;
            atmosph.fid = fopen([atmosph.expFileLocation,'\ATMOSPH.IN'],'r+');
            atmosph.C = textscan(atmosph.fid, '%s', 'delimiter', '','whitespace', '');
            atmosph.C = atmosph.C{1};
            fclose(atmosph.fid);
        end
       
        function setData(atmosph,param,paramValue,paramNum)
            for ii=1:length(atmosph.C)
                sLine = textscan(atmosph.C{ii}, '%s');
                ind = find(strcmp(param,sLine{1})==1);
                if ~isempty(ind) && ~strcmp(sLine{1}{1},'tAtm')
                    sLine = textscan(atmosph.C{ii+1}, '%s','whitespace','');
                    wline = writeLine(sLine{1},ind,num2str(paramValue));
                    atmosph.C{ii+1} = wline;
                    break
                elseif strcmp(sLine{1}{1},'tAtm')
                    numAtm = atmosph.countAtm(ii+1);
                    if numAtm > length(paramValue)
                        atmosph.C = [atmosph.C(1:ii+length(paramValue)); atmosph.C(end)];
                    elseif numAtm < length(paramValue)
                        sLine = textscan(atmosph.C{ii+1}, '%s','whitespace','');
                        newCells = repmat(sLine{1},[length(paramValue) 1]);
                        atmosph.C = [atmosph.C(1:ii); newCells ;atmosph.C(end)];
                    end
                    sLine = textscan(atmosph.C{ii}, '%s');
                    ind = find(strcmp(param,sLine{1})==1);
                    if ~isempty(ind)
                        kk = 0;
                        for jj=ii+1:ii+length(paramValue)
                            kk = kk + 1;
                            sLine = textscan(atmosph.C{jj}, '%s','whitespace','');
                            if length(ind) > 1
                                wline = writeLine(sLine{1},ind(paramNum),num2str(paramValue(kk)));
                            else
                                wline = writeLine(sLine{1},ind,num2str(paramValue(kk)));
                            end
                            atmosph.C{jj} = wline;
                        end
                        break
                    end
                end
            end
        end
                            
        function numAtm = countAtm(atmosph,ind)
            numAtm = 0;
            sLine = textscan(atmosph.C{ind+numAtm}, '%s');
            while ~(strcmp(sLine{1}{1},'end***') || strcmp(sLine{1}{1},'***'))
                numAtm = numAtm + 1;
                sLine = textscan(atmosph.C{ind+numAtm}, '%s');
            end
        end
      
        function update(atmosph)
            atmosph.fid = fopen([atmosph.expFileLocation,'\ATMOSPH.IN'],'w');
            frewind(atmosph.fid)
            for i=1:length(atmosph.C)
                wline = atmosph.C{i};
                fprintf(atmosph.fid,'%s\n',wline);
            end
            fclose(atmosph.fid);
            
            atmosph.fid = fopen([atmosph.expFileLocation,'\ATMOSPH.IN'],'r');
            atmosph.C = textscan(atmosph.fid, '%s', 'delimiter', '','whitespace', '');
            atmosph.C = atmosph.C{1};
            fclose(atmosph.fid);
            
        end % writeData   
    end % methods
end % classdef


        
        
        