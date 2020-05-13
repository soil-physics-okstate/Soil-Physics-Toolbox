Derek Groenendyk
4/23/2014
Notes on files for running HYDRUS in batch mode and other associated functions.


Update:

5/9/2012
I changed the results directory location in HYDRUS_Class, it saves the files
to a "Results" folder inside of your HYDRUS projects folder (as provided when 
initializing or creating the HYDRUS object).


Running a project:
( assumes you have already made a base case HYDRUS project using the GUI)
Required files - HYDRUS_Class.m, writeLine.m, readLine.m, IN_Class.m

Open HYDRUS_main.m
Provide: "exp" name, "exp" directory
Choose method: 1- single run, 2- multiple runs
Provide parameters headings and data if changing Selector.IN file information
Remember to set "numSims"
Run Matlab script


Notes:
I have files that modify the following files 
(not all written in Matlab however):

METEO.IN
ATMOSPH.IN (where the time-variable boundary conditions live)
SELECTOR.IN
PROFILE.OUT
BALANCE.OUT
T_LEVEL.OUT


