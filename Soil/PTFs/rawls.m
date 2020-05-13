function [theta_fc theta_wp] = rawls(Sand,Clay,OM)
%RAWLS is a simple linear pedotransfer function which estimates soil
%volumetric water content at field capacity and permanent wilting point
%based on sand, clay, and organic matter content. The reference is 
%
%Rawls, W.J., L.R. Ahuja, and D.L. Brakensiek.1992. Estimating soil 
%hydraulic properties from soils data. p. 329-340. In Proc. International 
%Workshop on Indirect Methods for Estimating the Hydraulic Properties of 
%Unsaturated Soils, Riverside, California. Univ. of California, Riverside.
%
%soil = [Sand Clay OM] all expressed in percentages
%theta_fc and theta_wp are expressed as volume fractions

theta_fc = 0.2576 -0.0020*Sand + 0.0036*Clay + 0.0299*OM; %-33 kPa, Table 1, Rawls et al. 1989. 
theta_wp = 0.0260 + 0.005*Clay +0.0158*OM; %-1500kPa, Table 1, Rawls et al. 1989
