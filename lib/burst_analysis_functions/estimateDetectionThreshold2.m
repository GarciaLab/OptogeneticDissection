% Script to scan through parameter space for selected biological variables
function [sweepInfo, param_val_vec] = estimateDetectionThreshold2(sweepInfo, param_val_vec, gillespie, use_flags)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% look at the lower piece
fluo_thresh = 1e5;
fluo_vec_raw = sweepInfo.raw_fluo_vec;
exp_fluo_vec = fluo_vec_raw(fluo_vec_raw<=fluo_thresh);
pd_fluo_vec = gillespie.fluo_ms2_array_noise(use_flags);%
pd_fluo_vec = pd_fluo_vec(pd_fluo_vec<=fluo_thresh);

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
options = optimoptions('lsqnonlin','Display','off');

% extract bounds and perform fits
paramList = sweepInfo.paramList;
ub_vec = sweepInfo.param_vals(end,strcmp(paramList,'HC_detect')|strcmp(paramList,'KD_detect'));
lb_vec = sweepInfo.param_vals(1,strcmp(paramList,'HC_detect')|strcmp(paramList,'KD_detect'));
guess_vec = (ub_vec + lb_vec)/2;
try
    fit = lsqnonlin(ob_fun,guess_vec,lb_vec,ub_vec,options);
catch
    fit(1) = sweepInfo.HC_detect;
    fit(2) = sweepInfo.KD_detect;    
end    

% generate predicted miss probability curve
% detection_stats.exp_detection_counts = exp_counts;
% detection_stats.exp_detection_counts_renorm = exp_counts_norm;
% detection_stats.predicted_detection_counts = pd_counts;
% detection_stats.detection_fluo_axis = fluo_axis_long;
sweepInfo.p_missed_detection_raw = p_missed_detection_long;
% sweepInfo.fluo_ref_curve = linspace(0,2*nanmax(fluo_vec_raw),1e3);
sweepInfo.p_miss_ref_vec = fit(2)^fit(1) ./ (fit(2)^fit(1) + sweepInfo.fluo_ref_curve.^fit(1));
param_val_vec(strcmp(paramList,'HC_detect')) = fit(1);
param_val_vec(strcmp(paramList,'KD_detect')) = fit(2);
