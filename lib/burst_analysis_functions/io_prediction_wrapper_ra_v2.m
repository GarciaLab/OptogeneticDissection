function sweepResults = io_prediction_wrapper_ra_v2(sweepInfo,sweepResults)
         
    param_val_vec = sweepResults.param_val_vec;    
    [sweepResults.param_val_vec, fluo_array_zeros, ra_cdf_pd_mean, ~, ~, gillespie]...
                                      = io_prediction_ra(param_val_vec, sweepInfo); 
    
    % calculate likelihood score
    delta_cdf = ra_cdf_pd_mean-sweepInfo.reactivation_cdf;
    ra_cdf_ste = sweepInfo.reactivation_cdf_ste;
    logL_ra_cdf = -(delta_cdf./ra_cdf_ste./2).^2 ;%- log(2*pi*ra_cdf_ste.^2);
    sweepResults.ra_fit = mean(logL_ra_cdf);
    sweepResults.ra_r2 = -nanmean((delta_cdf./nanmean(sweepInfo.reactivation_cdf)).^2);
    sweepResults.ra_N = sum(~isnan(delta_cdf));
    
    if sweepInfo.keep_prediction_flag
        fluo_out = fluo_array_zeros;
        for i = 1:size(fluo_out,2)
            start_i = find(fluo_out(:,i)>0,1);
            stop_i = find(fluo_out(:,i)>0,1,'last');
            fluo_out(1:start_i-1,i) = NaN;
            fluo_out(stop_i+1:end,i) = NaN;
        end

        % save basic info
        sweepResults.ra_time_cdf_predicted = ra_cdf_pd_mean;
        sweepResults.ms2_traces_observed_ra = fluo_out;
%         sweepResults.ms2_traces_true_ra = fluo_rarray;
        
        % save trace details        
        sweepResults.promoter_state_array_ra = gillespie.promoter_state_array;
        sweepResults.time_coarse_ra = gillespie.t_ref;
        sweepResults.time_fine_ra = gillespie.t_ref_full;                
        sweepResults.knirps_traces_ra = permute(gillespie.tf_ref_in,[1 3 2]);
%         sweepResults.reactivation_time_vec = reactivation_time_vec;
        sweepResults.kon_curves_ra = permute(gillespie.kon_curve_in,[1 3 2]); 
        sweepResults.koff_curves_ra = permute(gillespie.koff_curve_in,[1 3 2]); 

    end