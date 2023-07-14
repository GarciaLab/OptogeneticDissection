function [sweepResults, sweepInfo] = io_prediction_wrapper_cm(sweepInfo,sweepResults)
    
    % run through options 
    param_val_vec = sweepResults.param_val_vec;            
    
    [sweepResults.param_val_vec, sweepInfo, fluo_array_zeros, pd_fluo_curve, traceStruct] = io_prediction_cm(param_val_vec,sweepInfo);
    
    % calculate simple log Likelihood, assuming gaussian distributed errors
    fluo_diffs = pd_fluo_curve'-sweepInfo.mean_fluo_trend_cm;
    logL_vec = -fluo_diffs.^2 ./ (2*sweepInfo.ste_fluo_trend_cm.^2);% - log(2*pi*sweepInfo.ste_fluo_trend_cm.^2);
    sweepResults.fluo_io_fit = nanmean(logL_vec);  
    sweepResults.fluo_r2 = -nanmean((fluo_diffs./nanmean(sweepInfo.mean_fluo_trend_cm)).^2);  
    sweepResults.fluo_N = sum(~isnan(fluo_diffs));
    
    % save extra simulation info if option is flagged
    if sweepInfo.keep_prediction_flag
      
        fluo_out = fluo_array_zeros;
        
        sweepResults.fluo_io_curve_predicted = pd_fluo_curve';   
        sweepResults.ms2_traces_observed_cm = fluo_out;
        %sweepResults.ms2_traces_true_cm = fluo_array;
                
        sweepResults.promoter_state_array_cm = traceStruct.promoter_state_array;
        sweepResults.time_coarse_cm = traceStruct.t_ref;
        sweepResults.time_fine_cm = traceStruct.t_ref_full;        
        sweepResults.knirps_traces_cm = permute(traceStruct.tf_ref_in,[1 3 2]);          
        sweepResults.kon_curves_cm = permute(traceStruct.kon_curve_in,[1 3 2]);   
        sweepResults.koff_curves_cm = permute(traceStruct.koff_curve_in,[1 3 2]);   
       
    end