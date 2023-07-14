% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

% point to data directory
dataPath = '..\data\burst_analysis_data\';
DateStr = '03-Mar-2022 17-41-41';
outPath = [dataPath DateStr filesep];

% inisitalize sweep structure and set some parameters
sweepInfo.nParamIncrement = 10;
sweepInfo.keep_prediction_flag = false;
sweepInfo.r2_sweep_flag = false;
sweepInfo.koff_sweep_flag = true;

% set list of parameters to sample   
sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_detect','KD_detect', 'koff'                 ,'r2'};
sweepInfo.simFlags = [true,   true,    true,   false,       false,     sweepInfo.koff_sweep_flag, false];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sweep info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call wrapper script;
sweepInfo = initialize_io_sweep(dataPath, sweepInfo, DateStr);
                                    
% Specify range of parameter values to sweep
sweepInfo.param_vals = NaN(sweepInfo.nParamIncrement, length(sweepInfo.paramList));
sweepInfo.param_vals(:,1) = -linspace(4, 10, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,2) = linspace(3, 8, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,3) = linspace(sweepInfo.kon_75, 1.1*sweepInfo.kon_max, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,4) = linspace(0.5*sweepInfo.HC_detect, 2*sweepInfo.HC_detect, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,5) = linspace(0.5*sweepInfo.KD_detect, 2*sweepInfo.KD_detect, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,6) = linspace(sweepInfo.koff_25, sweepInfo.koff_75, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,7) = linspace(sweepInfo.r2_25, sweepInfo.r2_75, sweepInfo.nParamIncrement);

% add ofhter key parameters
sweepInfo.thresh_flags = contains(sweepInfo.paramList,'detect');

% sweepInfo.nFit = length(sweepInfo.paramList);
sweepInfo.nSim = sum(sweepInfo.simFlags);
sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nSim;

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);
sweepResults = initializeSweepValues(sweepInfo, sweepResults, []);              
    
tic
sweepResults = sweep_par_loop_v3(sweepInfo,sweepResults);    
toc

save([outPath 'sweepResults.mat'],'sweepResults')
save([outPath 'sweepInfo.mat'],'sweepInfo')

% Use MCMC to estimate parameter uncertainties
mcmcOptions = struct;
mcmcOptions.n_chains = 24; % number of unique chains to run
mcmcOptions.n_mcmc_steps = 2e3; % number of sampling steps
mcmcOptions.n_mcmc_prop_size = 0.04; % sets scale of proposal jump size
mcmcOptions.n_init = 100; % top N values from which to select starting point
mcmcOptions.n_traces = 50;
mcmcOptions.burn_in = 100;

tic
mcmcResults = io_mcmc_sampling(sweepInfo,sweepResults, mcmcOptions);
toc


save([outPath 'mcmcResults.mat'],'mcmcResults')
save([outPath 'mcmcOptions.mat'],'mcmcOptions')

