function [param_val_vec, sweepInfo, fluo_array_zeros, pd_fluo_curve, traceStruct] = io_prediction_cm(param_val_vec,sweepInfo)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = param_val_vec(strcmp(paramList,'HC_kon'));
    sweepInfo.KD = param_val_vec(strcmp(paramList,'KD_kon'));    
    sweepInfo.kon = param_val_vec(strcmp(paramList,'kon0'));
    sweepInfo.koff = param_val_vec(strcmp(paramList,'koff0'));
    sweepInfo.HC_koff = param_val_vec(strcmp(paramList,'HC_koff'));
    sweepInfo.KD_koff = param_val_vec(strcmp(paramList,'KD_koff'));
    sweepInfo.r2 = param_val_vec(strcmp(paramList,'r2'));

    % final model-building step
    sweepInfo = generate_full_model(sweepInfo);

    % extract info
    if isfield(sweepInfo,'n_traces_ra') && isfield(sweepInfo,'n_traces_cm')
        n_traces = sweepInfo.n_traces_cm;
    else
        n_traces = sweepInfo.n_traces;
    end

    % randomly draw tf profiles    
    s_indices = 1:size(sweepInfo.knirps_array_cm,2);%randsample(1:size(sweepInfo.knirps_array_cm,2),n_traces,true);
    tf_profile_array = permute(sweepInfo.knirps_array_cm(:,s_indices),[1 3 2]);

    %% call simulation function
    traceStruct = synthetic_rate_gillespie_io_v5(sweepInfo,tf_profile_array);    

    % estimate detection threshold implied by these parameter values
    use_flags = sweepInfo.ref_flag_array(:,s_indices);        
    [sweepInfo, param_val_vec] = estimateDetectionThreshold2(sweepInfo, param_val_vec, traceStruct, use_flags);    

    fluo_array = traceStruct.fluo_ms2_array;

    % apply threshold
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<0) = 0;

    % use reference curve to determine which time points are detected and
    % which are not
    df = sweepInfo.fluo_ref_curve(2) - sweepInfo.fluo_ref_curve(1);
    ref_vec = [sweepInfo.fluo_ref_curve-df/2 sweepInfo.fluo_ref_curve(end)+df];
    fluo_groups = discretize(fluo_array_zeros,ref_vec);
    miss_probs = sweepInfo.p_miss_ref_vec(fluo_groups);
    miss_status = miss_probs > rand(size(miss_probs));
    fluo_array_zeros(miss_status) = 0;    

    % apply time fiter to remove transient effects and enforce consistency
    % with experimental data
    time_filter = sweepInfo.time_axis_cm >= sweepInfo.time_range_cm(1) & ...
                  sweepInfo.time_axis_cm <= sweepInfo.time_range_cm(2);

    knirps_array_rs = permute(tf_profile_array,[1 3 2]);
    knirps_array_calc = knirps_array_rs(time_filter,:);
    fluo_array_calc = fluo_array_zeros(time_filter,:);

    % calculate predicted io function
    knirps_axis = sweepInfo.knirps_axis_cm;
    knirps_bins = linspace(2,12,length(knirps_axis)+1);
    pd_fluo_curve = NaN(size(sweepInfo.knirps_axis_cm));
    knirps_groups = discretize(permute(knirps_array_calc,[1 3 2]),knirps_bins);
    for k = 1:length(knirps_axis)
        pd_fluo_curve(k) = nanmean(fluo_array_calc(knirps_groups==k));
    end