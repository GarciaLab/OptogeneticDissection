function sweepInfo = initialize_io_sweep(resultsRoot, sweepInfo, cpHMMDataStr, varargin)

% set default values                                    
granularity = 1; % time resolution in seconds
n_traces = 5e2; % number of traces to simulate
max_time = 7*60; % max time to consider for reactivation cdf

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
sweepInfo.granularity = granularity;
sweepInfo.n_traces = n_traces;
sweepInfo.max_ra_time = max_time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% load markov system parameter info
sweepInfo = getMarkovSystemInfo_v2(sweepInfo,[resultsRoot cpHMMDataStr filesep]);

% record "true" profile
myCluster = parcluster('local');
max_workers = myCluster.NumWorkers;
sweepInfo.NumWorkers = min([24 max_workers]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate ground truth reference curves       
sweepInfo = addGroundTruthFields(sweepInfo, io_ref_ra, io_ref_cm, detection_stats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add parameter bounds
% Specify range of parameter values to sweep
sweepInfo.param_vals = NaN(sweepInfo.nParamIncrement, length(sweepInfo.paramList));

% set range of kon parameter values to explore
sweepInfo.param_vals(:,1) = -linspace(3, 8, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,2) = linspace(2, 6, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,3) = linspace(sweepInfo.kon0_lower, sweepInfo.kon0_upper, sweepInfo.nParamIncrement);

% set detection parameter values (these are only used to bound the fits)
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'HC_detect')) = linspace(0.5*sweepInfo.HC_detect, 2*sweepInfo.HC_detect, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'KD_detect')) = linspace(0.5*sweepInfo.KD_detect, 2*sweepInfo.KD_detect, sweepInfo.nParamIncrement);

% set koff bounds
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'koff0')) = linspace(0.07, 0.11, sweepInfo.nParamIncrement);%linspace(sweepInfo.koff0_lower, sweepInfo.koff0_upper, sweepInfo.nParamIncrement);
% sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'koffSlope')) = linspace(sweepInfo.koff_slope_lower, sweepInfo.koff_slope_upper, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'HC_koff')) = -linspace(0,4, sweepInfo.nParamIncrement);
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'KD_koff')) = linspace(2,6, sweepInfo.nParamIncrement);

% set r2 bounds
sweepInfo.param_vals(:,strcmp(sweepInfo.paramList,'r2')) = linspace(sweepInfo.r2_lower, sweepInfo.r2_upper, sweepInfo.nParamIncrement);

% add other key parameters
sweepInfo.thresh_flags = contains(sweepInfo.paramList,'detect');