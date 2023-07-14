% Script to test repeatability of calculated loss values
clear
close all

addpath(genpath('./lib'))

% point to data directory
dataPath = ['data' filesep  'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';

% outPath = [dataPath DateStr suffix filesep];
% mkdir(outPath)
figPath = ['.\fig\burst_analysis_results\consitency_checks' filesep];
mkdir(figPath)

% inisitalize sweep structure and set some parameters
sweepInfo.nParamIncrement = 5;
sweepInfo.keep_prediction_flag = false;
sweepInfo.r2_sweep_flag = false;
sweepInfo.koff_sweep_flag = true;
sweepInfo.koff_slope_sweep_flag = true;

% set list of parameters to sample   
sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_koff','koff0','r2','HC_detect','KD_detect'};
sweepInfo.simFlags = [true,   true,    true,    false(1,5)];
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'koff0')) = sweepInfo.koff_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'HC_koff')) = sweepInfo.koff_slope_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'r2')) = sweepInfo.r2_sweep_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sweep info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_traces_vec = [50 100 200 500 2e3];
master_struct = struct;

for n = 1:length(n_traces_vec)
  
    % call wrapper script;
    sweepInfo = initialize_io_sweep(dataPath, sweepInfo, DateStr,'n_traces',n_traces_vec(n));                                    

    % sweepInfo.nFit = length(sweepInfo.paramList);
    sweepInfo.nSim = sum(sweepInfo.simFlags);
    sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nSim;

    % initialize vectors to store results
    sweepResults = struct;
    [sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);
    sweepResults = initializeSweepValues(sweepInfo, sweepResults, []);              

    tic
    master_struct(n).sweepResults1 = sweep_par_loop_v3(sweepInfo,sweepResults);    
    toc

    tic
    master_struct(n).sweepResults2 = sweep_par_loop_v3(sweepInfo,sweepResults);    
    toc
end
  disp('Done.')
%% Test repeatability
close all
for n = 1:length(n_traces_vec)
    sweepResults1 = master_struct(n).sweepResults1;
    sweepResults2 = master_struct(n).sweepResults2;
    ra_loss_1 = -[sweepResults1.ra_r2];
    ra_loss_2 = -[sweepResults2.ra_r2];

    cm_loss_1 = -[sweepResults1.fluo_r2];
    cm_loss_2 = -[sweepResults2.fluo_r2];

    total_loss_1 = ra_loss_1 ./ -nanmean(ra_loss_1) + cm_loss_1 ./ -nanmean(cm_loss_1);
    total_loss_2 = ra_loss_2 ./ -nanmean(ra_loss_2) + cm_loss_2 ./ -nanmean(cm_loss_2);


    tot_loss = figure;
    scatter(total_loss_1,total_loss_2)

    xlabel('total loss (run 1)')
    ylabel('total loss (run 2)')
    grid on
    saveas(tot_loss,[figPath 'total_loss_scatter.png'])
    saveas(tot_loss,[figPath 'total_loss_scatter.pdf'])

    mdl = fitlm(total_loss_1,total_loss_2);
    master_struct(n).mdl = mdl;
end

