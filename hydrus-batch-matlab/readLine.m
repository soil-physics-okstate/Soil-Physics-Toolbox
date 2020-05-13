function [ paramStart,params,spaces ] = readLine( string )

    spaceLen = 0;
    paramLen = 0;
    paramStart = [];
    spaces = [];
    params = [];
    
    string = string{1};
    
    for ii=1:length(string)
       if ii == length(string)
           if strcmp(string(ii),' ')
               if ~strcmp(string(ii-1),' ')
                    params(end+1) = paramLen;
                    spaces(end+1) = spaceLen;
               end
               spaces(end+1) = (length(string)-(paramStart(end)+params(end) - 1));
           else
               if paramLen > 0
                   paramLen = paramLen + 1;
                   params(end+1) = paramLen;
                   spaces(end+1) = spaceLen;  %might not need, DGG 1/29/2015
               else
%                    spaces(end+1) = (length(string) - (paramStart(end) + params(end)) - 1);
                   spaces(end+1) = spaceLen;
                   paramStart(end+1) = ii;
                   params(end+1) = 1;
               end
           end
       else
           if strcmp(string(ii),' ')
               spaceLen = spaceLen + 1;
           else
               paramLen = paramLen + 1;
               if ii == 1
                   paramStart(end+1) = ii;
%                elseif ~strcmp(string(ii-1),' ')
%                    paramStart(end+1) = ii;
               elseif strcmp(string(ii-1),' ')
                   paramStart(end+1) = ii;
               end
               if ii+1 < length(string) - 1
                   if strcmp(string(ii+1),' ')
                       spaces(end+1) = spaceLen;
                       spaceLen = 0;
                       params(end+1) = paramLen;
                       paramLen = 0;
                   end
               end
           end
       end
    end
end
