function [model,param]=semivarfit3(d,V,N)
% function [model,param]=semivarfit3(d,V,N,map_date_str,depth)
% Removed map_date_str and depth jb to 4-27-2020
% SEMIVARFIT3 fits and compares variogram models to empirical variograms
% Inputs
%       d is the mean distance of each lag distance class
%       V is the vector of semivariances corresponding to the lag distances
%       N is the number of elements in each bin
%
% Outputs
%       model is a string of the chosen model
%       param is the corresponding parameters for the chosen model, which
%             is in the form of {nugget, sill, range}. 
%
% This function is adapted from semivarfit.m created by Andres Patrignani
% Modified by Jingnuo Dong 7/4/2015
% Modified by Tyson Ochsner 2/28/2020
% Removed gaussV model 2020-03-16

modelnamesgeo = {'spherV','exponV'};

% initial guess values for the parameters
param0 = [min(V), range(V)/2, mean(d)]; %initial guesses for nugget, sill, and range parameters
paramlist = cell(length(modelnamesgeo),1);
MSE = nan(length(modelnamesgeo),1);
lb = [0 0 0]; % lower bounds for nugget, sill, and range
ub = [max(V) max(V) 10*max(d)]; % upper bounds for nugget, sill, and range
options = optimoptions(@lsqnonlin,'Display','off'); % suppress command line outputs

% loop to fit each model in the list
for i=1:length(modelnamesgeo)
    fh = @(b)(feval(modelnamesgeo{i},b,d)-V).*N; % anonymous function computing the difference between modeled and measured variograms weighted by number of pairs in each bin
    [paramlist{i},resnorm] = lsqnonlin(fh,param0,lb,ub,options); % optimize semivariogram parameters, resnorm = weighted SSE
    MSE(i) = resnorm/length(d); % mean weighted squared error
end

% Chose the model with the minimum weighted MSE
[~,model_idx] = min(MSE);
param = paramlist{model_idx};
model = modelnamesgeo{model_idx};

% plot empirical variograms and the fitted models
figure
for i=1:length(modelnamesgeo)
    subplot(2,2,i)
    set(gca,'FontSize',14);
    plot(d,V,'ok')
    hold on
    k_alt = feval(modelnamesgeo{i},paramlist{i},d); %modeled semivariance
    plot(d,k_alt,'--b') %add model to plot
    ylabel('Semivariance ((cm^{-3} cm^{ 3})^{ 2})');
    xlabel('Lag distance (m)');
    title(modelnamesgeo{i})
    hold off
end
% print(strcat('../output/semivariogram/plots/semivariogram_',depth,'cm_', map_date_str), '-dpng');
% 
% % save model data
% dirOut = '../output/semivariogram/model/';
% fileName = strcat(dirOut, 'model_', depth, 'cm_', map_date_str, '.csv');
% fileOut = fopen(fileName, 'w');
% fprintf(fileOut, '%s\n', model);
% fclose(fileOut);
% 
% dlmwrite(fileName, param, 'precision', '%g', '-append');    

end
