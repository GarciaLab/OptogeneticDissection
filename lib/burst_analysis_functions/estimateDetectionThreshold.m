% Script to scan through parameter space for selected biological variables
function detection_stats = estimateDetectionThreshold(projectName, dataRoot, inferenceRoot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
resultsRoot = [dataRoot projectName filesep];
load([resultsRoot 'spot_struct.mat'],'spot_struct')

% get markov system parameters 
sweepInfo = struct;
sweepInfo = getMarkovSystemInfo_basic(sweepInfo,inferenceRoot);

% set basic parameters
sweepInfo.granularity = 1;
sweepInfo.n_traces = 100;
sweepInfo.seq_length = 120;
sweepInfo.deltaT = 20;
sweepInfo.granularity = 1; % nothing faster than a second
sweepInfo.keep_prediction_flag = true;
sweepInfo.tf_dependent_flags = false(2);
sweepInfo.RateMatrix = sweepInfo.R2_orig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sweepInfo.HC = 1;
sweepInfo.KD = 1;

tf_profile_array = rand(sweepInfo.seq_length,1,sweepInfo.n_traces);
    
gillespie = synthetic_rate_gillespie_io_v3(sweepInfo,tf_profile_array);        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% look at the lower piece
fluo_thresh = 1e5;
fluo_vec_raw = [spot_struct.fluo];
exp_fluo_vec = fluo_vec_raw(fluo_vec_raw<=fluo_thresh);
pd_fluo_vec = gillespie.fluo_ms2_array_noise(gillespie.fluo_ms2_array_noise<=fluo_thresh);


% figure(2);
fluo_bins = linspace(0,fluo_thresh,26);
fluo_axis = fluo_bins(1:end-1) + diff(fluo_bins)/2;

%% let's try to estimate a curve indicating probability of a missed
% detection as a function of true spot fluorescence
% should be 100% when F=0 and 0% (more or less) when F=1e5

exp_counts = histcounts(exp_fluo_vec,fluo_bins,'Normalization','Probability');
pd_counts = histcounts(pd_fluo_vec,fluo_bins,'Normalization','Probability');

% readjust probabilities
exp_counts_norm = exp_counts .* sum(pd_counts(end-2:end)) / sum(exp_counts(end-2:end)); % the idea is that, by the high end of the curve, missed counts should be negligible for exp 
p_missed_detection = (pd_counts - exp_counts_norm)./pd_counts;

% extend to reinforce need to congerge to zero
p_missed_detection_long = [p_missed_detection zeros(size(p_missed_detection))];
fluo_axis_long = [fluo_axis cumsum(diff(fluo_axis))+fluo_axis(end)];

% let's fit a simple hill curve to this
h_fun = @(x) x(2)^x(1) ./ (x(2)^x(1) + fluo_axis_long.^x(1));
ob_fun = @(x) p_missed_detection_long(1:end-1) - h_fun(x);
fit = lsqnonlin(ob_fun,[1 1e4],[-Inf -Inf],[Inf Inf]);

% generate predicted miss probability curve
detection_stats.exp_detection_counts = exp_counts;
detection_stats.exp_detection_counts_renorm = exp_counts_norm;
detection_stats.predicted_detection_counts = pd_counts;
detection_stats.detection_fluo_axis = fluo_axis_long;
detection_stats.p_missed_detection_raw = p_missed_detection_long;
detection_stats.fluo_ref_curve = linspace(0,2*nanmax(fluo_vec_raw),1e3);
detection_stats.KD_detect = fit(2);
detection_stats.HC_detect = fit(1);
detection_stats.p_miss_ref_vec = fit(2)^fit(1) ./ (fit(2)^fit(1) + detection_stats.fluo_ref_curve.^fit(1));

