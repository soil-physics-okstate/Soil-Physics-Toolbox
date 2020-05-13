% SELECTORIN.m
% Created by Derek Groenendyk
% 5/8/2012
% This a class that accesses the SELECTOR.IN file

classdef SELECTORIN < handle
    properties
        expFileLocation
        fid
        C
    end
      
    methods        
        function selector = SELECTORIN(directory)
            selector.expFileLocation = directory;
            selector.fid = fopen([selector.expFileLocation,'\SELECTOR.IN'],'r+');
            selector.C = textscan(selector.fid, '%s', 'delimiter', '','whitespace', '');
            selector.C = selector.C{1};
            fclose(selector.fid);
        end

        function setData(selector,param,paramValue,matNum)
            if (nargin < 4)  ||  isempty(matNum)
                matNum = 1;
            end
            for i=1:length(selector.C)
                idx = i;
                sLine = textscan(selector.C{idx}, '%s');
                ind = find(strcmp(param,sLine{1})==1);
                
                if idx == 26 && matNum > 1
                    numMat = selector.countMat();
                    if numMat < matNum
                        for j = numMat + 27:matNum+26
                            selector.C = [selector.C(1:j-1); selector.C(27); selector.C(j:end)];
                        end                        
                    end
                    idx = matNum+25;
                end
                
                if strcmp(param,'TPrint(1),TPrint(2),...,TPrint(MPL)')
                    ind = find(strcmp(selector.C,'TPrint(1),TPrint(2),...,TPrint(MPL)'));
                    endLine = find(strcmp(selector.C,'BLOCK D:')==1);
                    if isempty(endLine)
%                         endLine = find(strcmp(selector.C,'*** END OF INPUT FILE ''SELECTOR.IN'' ************************************')==1);
                        selector.C = [selector.C(1:ind+1); selector.C(end)];
                    else
                        selector.C = [selector.C(1:ind+1); selector.C(endLine:end)];
                    end
                    jj = 0;
                    idx = 0;
                    kk = 0;
                    while jj < length(paramValue)
                        sLine = textscan(selector.C{ind+1+kk}, '%s','whitespace','');
                        jj = jj + 1;
                        idx = idx + 1;
                        wline = writeLine(sLine{1},idx,num2str(paramValue(jj)));
                        selector.C{ind+1+kk} = wline;
                        if mod(jj,6) == 0 && jj+1 ~= length(paramValue)
                            kk = kk + 1;
                            idx = 0;
                            sLine = textscan(selector.C{ind+1}, '%s','whitespace','');
                            selector.C = [selector.C(1:ind+kk); sLine{1}; selector.C(ind+kk+1:end)];
                        end
                    end
                    if mod(jj,6) ~= 0 
                        for ii=6:-1:idx+1
                            sLine = textscan(selector.C{ind+1+kk}, '%s','whitespace','');
                            wline = writeLine(sLine{1},ii,'');
                            selector.C{ind+1+kk} = wline;
                        end
                    end
                    break
                end 
                
                if ~isempty(ind) && ~strcmp(param,'TPrint(1),TPrint(2),...,TPrint(MPL)')
                    sLine = textscan(selector.C{idx+1}, '%s','whitespace','');
                    wline = writeLine(sLine{1},ind,num2str(paramValue));
                    selector.C{idx+1} = wline;
                    break
                end  
            end
        end
        
        function numMat = countMat(selector)
            numMat = 0;
            sLine = textscan(selector.C{27+numMat}, '%s');
            while ~strcmp(sLine{1}{1},'***')
                numMat = numMat + 1;
                sLine = textscan(selector.C{27+numMat}, '%s');
            end
        end
            
        function update(selector)
            selector.fid = fopen([selector.expFileLocation,'\SELECTOR.IN'],'w');
            frewind(selector.fid)
            for i=1:length(selector.C)
                wline = selector.C{i};
                fprintf(selector.fid,'%s\n',wline);
            end
            fclose(selector.fid);
            
            selector.fid = fopen([selector.expFileLocation,'\SELECTOR.IN'],'r');
            selector.C = textscan(selector.fid, '%s', 'delimiter', '','whitespace', '');
            selector.C = selector.C{1};
            fclose(selector.fid);
            
        end % writeData   
    end % methods
end % classdef


        
        
        