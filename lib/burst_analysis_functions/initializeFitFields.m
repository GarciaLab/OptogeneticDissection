function [sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults)    
  
  % track fits to observables of interest
  for i = 1:sweepInfo.nIterations
    
      sweepResults(i).ra_fit = NaN;%(sweepInfo.nIterations,1);
      sweepResults(i).ra_r2 = NaN;
      sweepResults(i).ra_N = NaN;
      
      sweepResults(i).fluo_io_fit = NaN;
      sweepResults(i).fluo_r2 = NaN;
      sweepResults(i).fluo_N = NaN;
           
      if sweepInfo.keep_prediction_flag
          sweepResults(i).ra_time_cdf_predicted = NaN;
          sweepResults(i).fluo_io_curve_predicted = NaN;          
          
          sweepResults(i).promoter_state_array_cm = NaN;          
          sweepResults(i).time_coarse_cm = NaN;
          sweepResults(i).time_fine_cm = NaN;
          
          sweepResults(i).ms2_traces_observed_cm = NaN;
          sweepResults(i).ms2_traces_true_cm = NaN;
          sweepResults(i).knirps_traces_cm = NaN;          
          sweepResults(i).kon_curves_cm = NaN;   
          sweepResults(i).koff_curves_cm = NaN;   
          
          sweepResults(i).promoter_state_array_ra = NaN;
          sweepResults(i).time_coarse_ra = NaN;
          sweepResults(i).time_fine_ra = NaN;
          sweepResults(i).ms2_traces_observed_ra = NaN;
          sweepResults(i).ms2_traces_true_ra = NaN;
          sweepResults(i).knirps_traces_ra = NaN;
          sweepResults(i).reactivation_time_vec = NaN;
          sweepResults(i).kon_curves_ra = NaN;  
          sweepResults(i).koff_curves_ra = NaN;  
      end
  end