% PCAP HYDRUS batch run sample 
% Ben Paras
% Utilizes Derek's MATLAB code to run HYDRUS 1D. Modified to also run 2D. 
% 1/25/2016


run2D = 1;

%% HYDRUS1D
% Set up was done in GUI before modifying the rest of the files 

% Specify diretory and other parameters
% expDir_1 = 'C:\Users\Public\Documents\PC-Progress\Hydrus-1D 4.xx\Examples\Direct\PCAP_sim_2';
% h1d_calc hates spaces in path names so I copied the model into a path
% that doesn't have spaces 
if run2D == 0;
    expDir_1 = 'C:\Users\bkpar_000\Google_Drive\Thesis_bkparas\PCAP_sim_2';
    obs_nod_dir = [expDir_1 '\Obs_Node.out'];
    HYDRUS_1 = HYDRUS_Class(expDir_1);
    noCMDWindow = 1;
    n_sims = 3;
    q_sum = cell(3,1);
    % Use ATMOSPH.IN object 
    atmo_in_1 = ATMOSPHIN(expDir_1);
    Prec = [24.96 9.67 5.37];



    for i = 1 : n_sims
        atmo_in_1.setData('Prec', Prec(i),0);
        atmo_in_1.update();
        HYDRUS_1.run_hydrus(noCMDWindow,1);
        obs_nod = importdata(obs_nod_dir, ' ', 11);
        q_sum{i} = cumsum(obs_nod.data(:,4));
        HYDRUS_1.outputResults('PCAP_loam',num2str(i));
    end
    
else
    %% HYDRUS 2D/3D
    % Specify directory and other parameters
    expDir_2 = 'C:\Users\bkpar_000\Google_Drive\Thesis_bkparas\H3D2_PCAP_sim_2';
    q_dir = [expDir_2 '\v_Mean.out'];
    HYDRUS_2 = HYDRUS_Class(expDir_2);
    noCMDWindow = 1;
    n_sims = 6;
    q_top = cell(n_sims,1);
    q_plate = cell(n_sims,1);
    h_top = cell(n_sims,1);
    % Modify ATMOSPH.IN object
    atmo_in_2 = ATMOSPHIN(expDir_2);
    % rt = var flux 1 in HYDRUS2D/3D
    rt = [-9.67 -5.37 -3.21 -2.02 -1.32 -0.9054];
    paramNum = 1; % HYDRUS 2D/3D params in ATMOSPH.IN have multiple ones of the same name
    atmo_in_2.setData('rt',rt(1),paramNum);
    atmo_in_2.update();
    HYDRUS_2.run_hydrus(noCMDWindow,2);
    HYDRUS_2.outputResults('PCAP_loam_2D',num2str(1));
    q = importdata(q_dir, ' ' , 12); % 12 is the number of lines to skip before reading
    q_top{1} = q.data(:,6);
    q_plate{1} = q.data(:,7);
end