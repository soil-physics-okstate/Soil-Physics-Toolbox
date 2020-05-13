% PROFILEDAT.m
% Created by Derek Groenendyk
% 1/27/2015
% This a class that accesses the ATMOSPH.IN file

classdef PROFILEDAT < handle
    properties
        expFileLocation
        fid
        C
    end
      
    methods        
        function profile = PROFILEDAT(directory)
            profile.expFileLocation = directory;
            profile.fid = fopen([profile.expFileLocation,'\PROFILE.DAT'],'r+');
            profile.C = textscan(profile.fid, '%s', 'delimiter', '','whitespace', '');
            profile.C = profile.C{1};
            fclose(profile.fid);
        end
       
        function setData(profile,param,paramValue)
            if strcmp(param,'numLayers')
                row = 5;
                column = 1;                
                sLine = textscan(profile.C{row}, '%s','whitespace','');
                wline = writeLine(sLine{1},column,num2str(paramValue));
                profile.C{row} = wline;
            elseif strcmp(param,'columnLength')
                row = 4;
                column = 2;
                sData = sprintf('%e', paramValue);
                sLine = textscan(sData,'%8.7f %2s  %03d');
                paramValue = sprintf('-%fe+%03d',sLine{1},sLine{3});
                sLine = textscan(profile.C{row}, '%s','whitespace','');
                wline = writeLine(sLine{1},column,num2str(paramValue));
                profile.C{row} = wline;
            elseif strcmp(param,'obsLoc')
                endLine = profile.C{end};
                if length(strtrim(endLine)) == 1
                    profile.C = [profile.C(1:end-1) ; sprintf('% 5d',length(paramValue)) ; sprintf('% 5d',paramValue)];
                else
                    profile.C{end-1} = sprintf('% 5d',length(paramValue));
                    profile.C{end} = sprintf('% 5d',paramValue);
                end
            elseif ischar(param)                
                tempEndLine1 = profile.C{end};
                sLine = textscan(profile.C{end-1}, '%s');
                extraLineFlag = false;
                if length(strtrim(tempEndLine1)) ~= 1
                    tempEndLine1 = profile.C{end-1};
                    tempEndLine2 = profile.C{end};
                    sLine = textscan(profile.C{end-2}, '%s');
                    extraLineFlag = true;
                end
                numLayers = sLine{1}{1};
                sLine = textscan(profile.C{5}, '%s');
                ind = find(strcmp(sLine{1},param)==1)-3;
                if str2double(numLayers) < length(paramValue)
%                     'removing lines'
                    profile.C = profile.C(1:5+length(paramValue));
                elseif str2double(numLayers) < length(paramValue)
%                     'adding lines'
                    sLine = textscan(profile.C{6}, '%s','whitespace','');
                    newCells = repmat(sLine{1},[length(paramValue) 1]);
                    profile.C = [profile.C(1:5); newCells];
                end
                if str2double(numLayers) ~= length(paramValue)
%                     'adding end line'
                    if extraLineFlag
                        profile.C{end+1} = tempEndLine1;
                        profile.C{end+1} = tempEndLine2;
                    else
                        profile.C{end+1} = tempEndLine1;
                    end
                end
                
                for jj=6:5+length(paramValue)
                    sLine = textscan(profile.C{jj}, '%s','whitespace','');
                    wline = writeLine(sLine{1},1,sprintf('%d',jj-5));
                    profile.C{jj} = wline;                
                end
                
                [paramValue,exps] = profile.readExp(paramValue);
                if ~isempty(ind)
                    kk = 0;
                    for jj=6:5+length(paramValue)
                        kk = kk + 1;
                        sLine = textscan(profile.C{jj}, '%s','whitespace','');
                        
                        if exps(jj-5) >= 0.
                            val = sprintf('%fe+%03d',paramValue(jj-5),exps(jj-5));
                        elseif exps(jj-5) < 0.
                            val = sprintf('%fe%04d',paramValue(jj-5),exps(jj-5));
                        end
                        
                        wline = writeLine(sLine{1},ind,val);
                        profile.C{jj} = wline;
                    end
                end
            end
        end
        
        function [A,B] = readExp(profile,data)
            sData = sprintf('%e,',data);
            str1 = sprintf('%%%d.%df',8,6);
            results = textscan(sData,strcat(str1,' %2s  %03d'),'delimiter',',');
            A = results{1};
            B = results{3};
        end
            
        function update(profile)
            profile.fid = fopen([profile.expFileLocation,'\PROFILE.DAT'],'w');
            frewind(profile.fid)
            for i=1:length(profile.C)
                wline = profile.C{i};
                fprintf(profile.fid,'%s\n',wline);
            end
            fclose(profile.fid);
            
            profile.fid = fopen([profile.expFileLocation,'\PROFILE.DAT'],'r');
            profile.C = textscan(profile.fid, '%s', 'delimiter', '','whitespace', '');
            profile.C = profile.C{1};
            fclose(profile.fid);
            
        end % writeData   
    end % methods
end % classdef


        
        
        