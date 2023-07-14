% Script to scan through parameter space for selected biological variables
function sweepInfo = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,...
                                      paramVals,keep_prediction_flag,varargin)

% set default values                                    
granularity = 1; % time resolution
n_traces = 100; % number of traces to simulate
max_time = 7*60;

for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

% load data
load([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
load([resultsRoot 'io_ref_cm.mat'],'io_ref_cm')
load([resultsRoot 'detection_stats.mat'],'detection_stats')

% set basic parameters
sweepInfo = struct;
sweepInfo.nParamIncrement = nParamIncrement;
sweepInfo.granularity = granularity;
sweepInfo.n_traces = n_traces;
sweepInfo.keep_prediction_flag = keep_prediction_flag;
sweepInfo.max_ra_time = max_time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% save simulation type
sweepInfo.simType = simType;

% load markov system parameter info
sweepInfo = getMarkovSystemInfo(sweepInfo,resultsRoot);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call Sweep function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

% record "true" profile
myCluster = parcluster('local');
max_workers = myCluster.NumWorkers;
sweepInfo.NumWorkers = min([24 max_workers]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate ground truth reference curves       
sweepInfo = addGroundTruthFields(sweepInfo, io_ref_ra, io_ref_cm, detection_stats);

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults, paramVals);             
sweepResults = initializeSweepValues(sweepInfo, sweepResults, paramVals);              

% call parallel sweep script
sweepInfo.NumWorkers = min([min([28 max_workers]), length(sweepResults)]);
tic
sweepResults= sweep_par_loop_v3(sweepInfo,sweepResults);    
toc

% recombine 
fnames = fieldnames(sweepResults);
if ~keep_prediction_flag
    for f = 1:length(fnames)
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
else
    for f = [1:11 length(fnames)]
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
    for f = 12:length(fnames)-1
        sweepInfo.(fnames{f}) = cat(3,sweepResults.(fnames{f}));
    end
end
clear sweepResults    
if isempty(paramVals)
    save([outPath 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
end