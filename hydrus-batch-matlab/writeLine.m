function [ wline ] = writeLine( string, index, param, lspaces )
%     param = num2str(param);
    flag = 0;
    if nargin > 3
        flag = 1;
        [starts, paramLens, whitespaces] = readLine(string);
        whitespaces(1) = whitespaces(1)-length(param)+1;
    else
        [starts, paramLens, whitespaces] = readLine(string);
    end
    starts
    paramLens
    whitespaces
    sLine = textscan(string{1}, '%s');
    sLine = sLine{1};
    oldLen = paramLens(index);
    newLen = length(param);
    
    if newLen > oldLen + whitespaces(index)
        param = sprintf('%e', str2num(param))
        newLen = length(param)
    end
    
    paramLens(index) = newLen;
    sLine(index) = cellstr(param);

    if index ~= length(sLine) && ~flag % double check
        whitespaces(index+1) = whitespaces(index+1) + oldLen - newLen;
    elseif ~flag
        if length(whitespaces) > length(paramLens)
            whitespaces(index+1) = whitespaces(index+1) + oldLen - newLen; 
        end
    end
    
    if starts(1) == 0
        wline = '';
        for i=2:length(whitespaces)
            wline = [wline , sLine{i-1},blanks(whitespaces(i))];
        end
        if length(whitespaces)-1 < length(paramLens)
            wline = [wline , sLine{end}];
        end
    else
        wline = '';
        for i=1:length(paramLens)
            wline = [wline , blanks(whitespaces(i)), sLine{i}];
        end
        if length(whitespaces) > length(paramLens)
            wline = [wline , blanks(whitespaces(end))];
        end
    end    
end

