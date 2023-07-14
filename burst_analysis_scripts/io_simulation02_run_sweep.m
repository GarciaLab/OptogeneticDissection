% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../lib'))

% point to data directory
dataPath = ['..' filesep 'data' filesep 'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_FINAL';
% suffix = '_RR2';
outPath = [dataPath DateStr suffix filesep];
mkdir(outPath)

% seed random number generator
rng(237);

% inisitalize sweep structure and set some parameters
sweepInfo.nParamIncrement = 15;
sweepInfo.keep_prediction_flag = false;
sweepInfo.r2_sweep_flag = false;
sweepInfo.kon_sweep_flag = false; 
sweepInfo.koff_sweep_flag = false;
if contains(suffix,'kon_only')
    sweepInfo.koff_kd_sweep_flag = false;
    sweepInfo.koff_slope_sweep_flag = false;
else
    sweepInfo.koff_kd_sweep_flag = true;
    sweepInfo.koff_slope_sweep_flag = true;
end
% set list of parameters to sample   
sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_koff','KD_koff','koff0','r2','HC_detect','KD_detect'};
sweepInfo.simFlags = [true,   true,   false(1,7)];
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'koff0')) = sweepInfo.koff_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'HC_koff')) = sweepInfo.koff_slope_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'KD_koff')) = sweepInfo.koff_kd_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'kon0')) = sweepInfo.kon_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'r2')) = sweepInfo.r2_sweep_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sweep info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call wrapper script;
sweepInfo = initialize_io_sweep(dataPath, sweepInfo, DateStr);                                    

% sweepInfo.nFit = length(sweepInfo.paramList);
sweepInfo.nSim = sum(sweepInfo.simFlags);
sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nSim;
sweepInfo.kon_only_flag = false;
if ~sweepInfo.koff_slope_sweep_flag && ~sweepInfo.koff_kd_sweep_flag
    sweepInfo.kon_only_flag = true;
end

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);
sweepResults = initializeSweepValues(sweepInfo, sweepResults, []);              
   
tic
sweepResults = sweep_par_loop_v3(sweepInfo,sweepResults);    
toc

save([outPath 'sweepResults.mat'],'sweepResults')
save([outPath 'sweepInfo.mat'],'sweepInfo')
