function [param_val_vec, fluo_array_zeros, ra_cdf_pd_obs, ra_cdf_pd_true, kon_curve_ra, gillespie] ...
                                = io_prediction_ra(param_val_vec, sweepInfo)

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
        n_traces = sweepInfo.n_traces_ra;
    else
        n_traces = sweepInfo.n_traces;
    end                                                       
    
    % randomly draw tf profiles    
    s_indices = 1:size(sweepInfo.knirps_array_ra,2);%(srandsample(1:size(sweepInfo.knirps_array_ra,2),n_traces,true);
    tf_profile_array = permute(sweepInfo.knirps_array_ra(:,s_indices),[1 3 2]);

    %% call simulation function
    gillespie = synthetic_rate_gillespie_io_v5(sweepInfo,tf_profile_array);

    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;   
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<0) = 0;
    
    % use reference curve to determine which time points are detected and
    % which are not    
    df = sweepInfo.fluo_ref_curve(2) - sweepInfo.fluo_ref_curve(1);
    ref_vec = [sweepInfo.fluo_ref_curve-df/2 sweepInfo.fluo_ref_curve(end)+df];
    fluo_groups = discretize(fluo_array_zeros,ref_vec);
    miss_probs = sweepInfo.p_miss_ref_vec (fluo_groups);
    miss_status = miss_probs > rand(size(miss_probs));
    fluo_array_zeros(miss_status) = 0;       
    
    % calculate stats for fraction of traces that actually turn off
    perturbation_frame = ceil(size(fluo_array,1)/2);
    window_size = perturbation_frame-1;
    index_vec = -window_size:window_size;
    off_frames = sweepInfo.off_frame_ref;
    
    % was trace on before and after perturbation?
    active_indices = 1*(fluo_array_zeros>0) .* index_vec';
    first_i_vec = min(active_indices);
    max_time = max(sweepInfo.reactivation_time);
    before_flags = first_i_vec < 0;
    
    % was trace OFF immediately preceding ON perturbation?
    off_flags = all(fluo_array_zeros(ismember(index_vec,off_frames),:)==0,1);
    
    % get predicted reactivation curve
    ra_cdf_pd_obs = calculate_ra_cdf(off_flags, before_flags, fluo_array_zeros, index_vec, max_time, sweepInfo, 0);
    
    % get "true" reactivation curve
    f_baseline = min(min(fluo_array(10:end,:)));
    ra_cdf_pd_true = calculate_ra_cdf(off_flags, before_flags, fluo_array, index_vec, max_time, sweepInfo, f_baseline);
    
    % get predicted mean kon curve
    kon_curve_ra = nanmean(permute(gillespie.kon_curve_in,[1 3 2]),2); 
    