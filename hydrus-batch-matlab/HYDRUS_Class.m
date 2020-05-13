% HYDRUS_Class.m
% Derek Groenendyk
% Created 5/4/2012
% Class that calls and runs HYRDUS in Matlab
% modified 1/27/2016 by Ben Paras

classdef HYDRUS_Class
    properties
        expFileLocation;
        D;
        fileExtList = {'.OUT' '.IN' '.TXT' '.DAT'};
    end
      
    methods
        function hydrus = HYDRUS_Class(directory)
            hydrus.expFileLocation = directory;
            hydrus.D = 1;
        end
        function run_hydrus(hydrus,noCMDWindow,D)
            selectIN = SELECTORIN(hydrus.expFileLocation);
%             if nargin > 1
%                 
%             end
            % Fill in unset optional values.
            switch nargin
                case 1
                    hydrus.D = D;
            end

            if noCMDWindow == 1
                selectIN.setData('lEnter','f');
                selectIN.update()
                command = ['h' num2str(hydrus.D) 'd_calc.exe ' hydrus.expFileLocation '> NUL 2>&1'];
                system(command);                
            else
                selectIN.setData('lEnter','t');
                selectIN.update()
                command = ['h' num2str(hydrus.D) 'd_calc.exe ' hydrus.expFileLocation ];
                [status,result] = system(command,'-echo');
            end
        end
        
        function outputResults(hydrus,folder,trial)
            if nargin < 2
                trial =  '1';
            end
            root = pwd;
            cd(hydrus.expFileLocation)
            listing = dir(pwd);
            cd(root)
                
            for ii=1:length(listing)-2
                [pathstr, name, ext] = fileparts(listing(ii+2).name);
                temp = cell2mat(strfind(hydrus.fileExtList,upper(ext)));
                
                if length(temp==1) > 0
                    INFILES{ii} = listing(ii+2).name;
                end
            end
            ind = strfind(hydrus.expFileLocation,'\');
            resultsDir = [hydrus.expFileLocation(1:ind(end)) 'Results\' folder '\Trial= ' trial];
            
            if exist(resultsDir,'dir') == 0
                mkdir(resultsDir)
            end
               
%         Will automatically overwrite files
          % If INFILES has empty elements, clear (HYDRUS2D will have some
          % because of extra folders it has) 
          INFILES = INFILES(~cellfun('isempty',INFILES));
            for i=1:length(INFILES)
                srcDir = [hydrus.expFileLocation '\' char(INFILES(i))];
%                 untested: code to prevent overwrite
%                if exist(srcDir,'file') > 0 
%                    ['File ' char(INFILES(i)) ' Overwritten']
%                end
                desDir = [resultsDir];
                copyfile(srcDir,desDir,'f');
            end
        end %outputResults
        
    end %methods
end %HYDRUS_Class
