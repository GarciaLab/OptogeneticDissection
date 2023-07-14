function mcmcOptions = getMCMCPriors_v2(mcmcOptions, sweepInfo, param_array_resamp)

    % load inference results
%     load([cpHMMPath filesep 'inference_summary.mat'],'inference_summary')

    mcmcOptions.memory = 7;
    mcmcOptions.deltaT = 20;
    mcmcOptions.t_MS2 = 1.4;
    
    % initialize vectors to store prior mean and sigma info
    mcmcOptions.prior_mean_vec = NaN(size(sweepInfo.paramList));
    mcmcOptions.prior_std_vec = NaN(size(sweepInfo.paramList));
   
    % for values that were swept, take mean and std of results as prior
    mcmcOptions.prior_mean_vec(sweepInfo.simFlags) = nanmean(param_array_resamp(:,sweepInfo.simFlags),1);
    mcmcOptions.prior_std_vec(sweepInfo.simFlags) = nanstd(param_array_resamp(:,sweepInfo.simFlags),[],1);
    
    % otherwise take default values and sigmas (kon0 and koff0 reflect
    % cpHMM resuts)
    mcmcOptions.prior_mean_vec(~sweepInfo.simFlags) = sweepInfo.defaultValues(~sweepInfo.simFlags);
    mcmcOptions.prior_std_vec(~sweepInfo.simFlags) = sweepInfo.priorSigmas(~sweepInfo.simFlags);   
    
    mcmcOptions.prior_mean_vec = mcmcOptions.prior_mean_vec(~sweepInfo.thresh_flags);
    mcmcOptions.prior_std_vec = mcmcOptions.prior_std_vec(~sweepInfo.thresh_flags);
    
    % adjust parameter bounds
    paramBounds = mcmcOptions.paramBounds;
    paramBounds(1,1:7) = mcmcOptions.prior_mean_vec-3*mcmcOptions.prior_std_vec;
    paramBounds(2,1:7) = mcmcOptions.prior_mean_vec+3*mcmcOptions.prior_std_vec;
    
    neg_flags = paramBounds(1,:)<0;
    neg_flags([1 4]) = false;
    paramBounds(1,neg_flags) = 0;
    mcmcOptions.paramBounds = paramBounds;
    
    % select inital points from sweep array
    mcmcOptions.paramInit = param_array_resamp(randsample(1:length(param_array_resamp),mcmcOptions.n_chains,true),1:7);
    
    % deal with "kon only" option
    if mcmcOptions.kon_only_flag
        %%%%
        % reset parameter initializations
        % set koff Kd to nearly zero. This eliminates Knirps-dependence
        koff_kd_ind = 5;
        mcmcOptions.paramInit(:,koff_kd_ind) = 1e-10;
        % set H to 1 for simplicity (choice does not really matter)
        koff_h_ind = 4;
        mcmcOptions.paramInit(:,koff_h_ind) = 1;
        
        %%%% update parameter means and param bounds
        mcmcOptions.paramBounds(:,[koff_h_ind, koff_kd_ind]) = repmat([1 1e-10],2,1);
        mcmcOptions.prior_mean_vec([koff_h_ind, koff_kd_ind]) = [1 1e-10];
        mcmcOptions.prior_std_vec([koff_h_ind, koff_kd_ind]) = [1 1]; % this does not matter. LogL contribution will be identical across models and will drop out
    end
    