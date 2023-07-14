function sweepResults = initializeSweepValues(sweepInfo, sweepResults, param_fit_array_sim) 
        
  simFlags = sweepInfo.simFlags;   
  
  % initialize new array to store actual results  
  if isempty(param_fit_array_sim)
      param_fit_array_sim = NaN(sweepInfo.nIterations, length(sweepInfo.simFlags)); 
      param_fit_array_sim(:,~sweepInfo.simFlags) = repmat(sweepInfo.defaultValues(~sweepInfo.simFlags),sweepInfo.nIterations,1);     
      
      iter = 1;
      for i = find(simFlags)
          vec_orig = sweepInfo.param_vals(:,i);
          vec1 = repmat(vec_orig,sweepInfo.nIterations / sweepInfo.nParamIncrement^iter,1);
          vec2 = repelem(vec1,sweepInfo.nIterations / sweepInfo.nParamIncrement^(sweepInfo.nSim-iter+1));     
          param_fit_array_sim(:,i) = vec2;
          iter = iter + 1;
      end 
      
      if sweepInfo.kon_only_flag
          param_fit_array_sim(:,(strcmp(sweepInfo.paramList,'KD_koff'))) = 1e-10; % make this essentially 0
          param_fit_array_sim(:,(strcmp(sweepInfo.paramList,'HC_koff'))) = 1; % the precise ßvalue is of little consequence heree
      elseif ~sweepInfo.simFlags(strcmp(sweepInfo.paramList,'KD_koff'))
          param_fit_array_sim(:,(strcmp(sweepInfo.paramList,'KD_koff'))) = param_fit_array_sim(:,(strcmp(sweepInfo.paramList,'KD_kon')));
      end
  end
  
  % assign  to structure
  for j = 1:size(param_fit_array_sim, 1)
      sweepResults(j).param_val_vec = param_fit_array_sim(j,:);
  end  