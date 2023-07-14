function [sweepInfo, sweepResults] = initializeFitFields_v2(sweepInfo,sweepResults,param_fit_array)
  
%   sweepResults = struct;
  
  % set list of parameters to sample   
  sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_detect','KD_detect'};
%   sweepInfo.fitFlags = [1 1 1 1 1];
%   sweepInfo.trueVals = [7,6e5,0,0,sweepInfo.R2_orig(2,1),sweepInfo.R2_orig(1,2)]; 
%   
%   if contains(sweepInfo.simType,'2')      
%       ka_index = strcmp(sweepInfo.paramList,'ka');
%       ks_index = strcmp(sweepInfo.paramList,'ks');
%       sweepInfo.fitFlags(ks_index | ka_index) = 0; 
%       if contains(sweepInfo.simType,'koff')
%           sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'kon')) = 0;
%       elseif contains(sweepInfo.simType,'kon')
%           sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;
%       end
%   elseif strcmp(sweepInfo.simType,'match_exp')
%       sweepInfo.fitFlags(1:end) = false;      
%   else
%       % koff is always fixed for 3 state models
%       sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;          
%   end
%   sweepInfo.nFit = sum(sweepInfo.fitFlags);
%   if strcmp(sweepInfo.simType,'match_exp')
%       sweepInfo.nIterations = 1;
%       sweepInfo.nParamIncrement = 1;
%   elseif ~isfield(sweepInfo, 'nIterations')
%       sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;
%   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set bounds on parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set sweep sampling hyperparameters
  sweepInfo.param_values = NaN(2,sweepInfo.nParamIncrement);
  
  % kon hill coefficient
  sweepInfo.param_bounds(:,1) = linspace(4, 8, sweepInfo.nParamIncrement);
  
  % kon KD
  sweepInfo.param_bounds(:,2) = linspace(4, 12, sweepInfo.nParamIncrement); 
  
  % kon max
  sweepInfo.param_bounds(:,3) = linspace(0.5*sweepInfo.kon_max, 2*sweepInfo.kon_max, sweepInfo.nParamIncrement); 

  if ~isempty(param_fit_array)
      sweepInfo.nIterations = size(param_fit_array,1);
  end
  % track fits to observables of interest
  for i = 1:sweepInfo.nIterations      
      
      sweepResults(i).ra_fit = NaN;
      sweepResults(i).fluo_io_fit = NaN;
      
      if sweepInfo.keep_prediction_flag
          sweepResults(i).ra_time_cdf_predicted = NaN;%(1,length(sweepInfo.reactivation_time));          
          sweepResults(i).fluo_io_predicted = NaN;
                              
          sweepResults(i).ms2_traces_observed_cm = NaN;
          sweepResults(i).ms2_traces_true_cm = NaN;
          sweepResults(i).knirps_traces_cm = NaN;          
          sweepResults(i).kon_curve_cm = NaN;                              
          sweepResults(i).kon_curves_cm = NaN;           
          
          sweepResults(i).ms2_traces_observed_ra = NaN;
          sweepResults(i).ms2_traces_true_ra = NaN;
          sweepResults(i).knirps_traces_ra = NaN;
          sweepResults(i).reactivation_time_vec = NaN;
          sweepResults(i).kon_curves_ra = NaN;           
      end
  end