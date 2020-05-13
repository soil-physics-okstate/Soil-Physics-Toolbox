%% Tutorial Block kriging
%
% IN PROCESS...

% Known data. S[x y z]
S = [61 139 477;
     63 140 696;
     64 129 227;
     68 128 646;
     71 140 606;
     73 141 791;
     75 128 783];

% Block
B = [64 132;
     64 138;
     70 138;
     70 132];

% Points within block
P = [66 134; % Point A
     66 136;  % Point B
     68 134;  % Point C
     68 136]; % Point D
 
% Visualize data
figure
scatter([S(:,1);P(:,1)] , [S(:,2);P(:,2)], 'ok') % Plot
hold on % Keep figure active to receive more information
minx = min(B(:,1)); % X coordinate of the lower-left corner
miny = min(B(:,2)); % Y coordinate of the lower-left corner
dx = max(B(:,1))-min(B(:,1));
dy = max(B(:,2))-min(B(:,2));
rectangle('Position',[minx,miny,dx,dy])
 
% krging
obs = S;
pred = P;
beta = [10 3.33]; % Assume the following semivariogram model
maxpoints = 7; % All of them
maxdist = 50;
Gmodel = 'exponV'; % Exponential semivariogram model
Kmodel = 'ordinary'; % Kriging model
[Zpred,Vpred] = krig(obs,pred,beta,maxpoints,maxdist,Gmodel,Kmodel);

