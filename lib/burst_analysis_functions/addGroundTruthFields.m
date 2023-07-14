function sweepInfo = addGroundTruthFields(sweepInfo,io_ref_ra,io_ref_cm,detection_stats)

    % RA fields. Take only times up to designated point    
    use_flags = io_ref_ra.reactivation_time_axis<=sweepInfo.max_ra_time;
    sweepInfo.reactivation_cdf = io_ref_ra.reactivation_time_cdf(use_flags);
    sweepInfo.reactivation_cdf_ste = io_ref_ra.reactivation_time_cdf_ste(use_flags);
    sweepInfo.reactivation_time = io_ref_ra.reactivation_time_axis(use_flags);
    sweepInfo.off_frame_ref = io_ref_ra.off_frame_ref;   
    sweepInfo.knirps_array_ra = io_ref_ra.knirps_array;  
    sweepInfo.fluo_array_ra = io_ref_ra.fluo_array;  
    sweepInfo.time_axis_ra = io_ref_ra.time_vec';    
    
    % fluo trend 
    sweepInfo.time_axis_cm = io_ref_cm.time_axis;   
    sweepInfo.raw_fluo_vec = io_ref_cm.raw_fluo_vec;   
    sweepInfo.ref_flag_array = io_ref_cm.ref_flag_array_ft;
    sweepInfo.ste_fluo_trend_cm = io_ref_cm.ste_fluo_trend;
    rm_flags = isnan(sweepInfo.ste_fluo_trend_cm) | sweepInfo.ste_fluo_trend_cm<0.2e4;
    sweepInfo.ste_fluo_trend_cm = sweepInfo.ste_fluo_trend_cm(~rm_flags);
    sweepInfo.mean_fluo_trend_cm = io_ref_cm.mean_fluo_trend(~rm_flags);    
    sweepInfo.knirps_axis_cm = io_ref_cm.mean_knirps_trend(~rm_flags);
    sweepInfo.knirps_array_cm = io_ref_cm.knirps_array_ft; 
    sweepInfo.fluo_array_cm = io_ref_cm.fluo_array;  
    sweepInfo.time_range_cm = io_ref_cm.time_range;  
    
    % add detection threshold info
    fnames = fieldnames(detection_stats);
    for f = 1:length(fnames)
        sweepInfo.(fnames{f}) = detection_stats.(fnames{f});
    end
    