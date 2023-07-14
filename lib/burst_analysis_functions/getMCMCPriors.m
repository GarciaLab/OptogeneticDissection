function mcmcOptions = getMCMCPriors(mcmcOptions, sweepInfo, cpHMMPath,param_array_resamp)

    % load inference results
    load([cpHMMPath filesep 'inference_summary.mat'],'inference_summary')

    mcmcOptions.memory = 7;
    mcmcOptions.deltaT = 20;
    mcmcOptions.t_MS2 = 1.4;
    
    % initialize vectors to store prior mean and sigma info
    mcmcOptions.prior_mean_vec = NaN(size(sweepInfo.paramList));
    mcmcOptions.prior_std_vec = NaN(size(sweepInfo.paramList));
   
    % for values that were swept, take mean and std of results as prior
    mcmcOptions.prior_mean_vec(sweepInfo.simFlags) = nanmean(param_array_resamp(:,sweepInfo.simFlags),1);
    mcmcOptions.prior_std_vec(sweepInfo.simFlags) = nanstd(param_array_resamp(:,sweepInfo.simFlags),[],1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OFF RATE
    % when we're not sweeping koff slope, use mean and quartiles as bounds
    koff_vec = [inference_summary.koff_mean]/60;
    koff_ste_vec = [inference_summary.koff_ste]/60;
    
    if ~sweepInfo.koff_slope_sweep_flag                
        koff0_guess = nanmean(koff_vec);
        koff0_std = nanstd(koff_vec);
    else
        % otherwise, conduct simple bootstrapped linear fit to obtain
        % bounds on intercept and slope
        knirps_vec = [inference_summary.knirps_mean];
        knirps_ste_vec = [inference_summary.knirps_ste];
        knirps_axis = linspace(0,1.1*max(knirps_vec));
        nan_ft = ~isnan(koff_vec);
        koff_weight_vec = (koff_ste_vec + 0.01*nanmean(koff_vec)).^-2; % extra factor prevents over-fitting of low-variance points
        nBoots = 100;
        [koff_trend, ~, ~, ~, param_95, param_05, ~]...
                                              = fit_lin_trend(...
                                              knirps_axis,knirps_vec(nan_ft),knirps_ste_vec(nan_ft), koff_vec(nan_ft),...
                                                koff_ste_vec(nan_ft),koff_weight_vec(nan_ft).^2,nBoots);        
        koff0_upper = param_95(1);
        koff0_lower = param_05(1);
        koff0_guess = koff_trend(1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ON RATE%
    % when we're not sweeping koff slope, use mean and quartiles as bounds
    kon_vec = [inference_summary.kon_mean]/60;
    kon_ste_vec = [inference_summary.kon_ste]/60;
        
    % otherwise, conduct simple bootstrapped linear fit to obtain
    % bounds on intercept and slope        
    nan_ft = ~isnan(kon_vec);
    kon_weight_vec = (kon_ste_vec + 0.01*nanmean(kon_vec)).^-2; % extra factor prevents over-fitting of low-variance points    
    [kon_trend, ~, ~, ~, kon_95, kon_05, ~]...
                                          = fit_lin_trend(...
                                          knirps_axis,knirps_vec(nan_ft),knirps_ste_vec(nan_ft), kon_vec(nan_ft),...
                                            kon_ste_vec(nan_ft),kon_weight_vec(nan_ft).^2,nBoots);        
    kon0_upper = kon_95(1);
    kon0_lower = kon_05(1);
    kon0_guess = kon_trend(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initiation rate bounds
    r1_mean = nanmean([inference_summary.r1_mean]);
    r2_guess = nanmean([inference_summary.r2_mean]);
    r2_lower = prctile([inference_summary.r2_mean],25);
    r2_upper = prctile([inference_summary.r2_mean],75);
    
    % noise in MS2 signal
    noise_mean = nanmedian([inference_summary.noise_mean]);        
    
    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    sweepInfo.R2_orig = [-kon0_guess  koff0_guess; 
                          kon0_guess -koff0_guess];    
    
    sweepInfo.r_emission = [r1_mean r2_guess]*sweepInfo.deltaT;
    sweepInfo.pi0 = [koff0_guess/(kon0_guess +koff0_guess) kon0_guess /(kon0_guess +  koff0_guess)];    
    sweepInfo.noise = noise_mean;
    
    sweepInfo.kon0_lower = kon0_lower;
    sweepInfo.kon0_upper = kon0_upper;
    sweepInfo.kon_mean = kon0_guess;
    
    sweepInfo.koff0_lower = koff0_lower;
    sweepInfo.koff0_upper = koff0_upper;
    sweepInfo.koff_mean = koff0_guess;
    
    sweepInfo.r2_lower = r2_lower*sweepInfo.deltaT;
    sweepInfo.r2_mean = r2_guess*sweepInfo.deltaT;
    sweepInfo.r2_upper = r2_upper*sweepInfo.deltaT;
    
    % r2 (we can be pretty strict about this one)
    mcmcOptions.r2_prior_mean = nanmean([inference_summary(use_inds).r2_mean])*mcmcOptions.deltaT;
    mcmcOptions.r2_prior_std = nanstd([inference_summary(use_inds).r2_mean])*mcmcOptions.deltaT;

    % koff0 (should be higher than most if not all cpHMM results)    
    koff_guess = prctile(koff_high_vec,95);
    mcmcOptions.koff0_prior_mean = koff_guess;
    mcmcOptions.koff0_prior_std = koff_guess - nanmean(koff_vec);%nanstd([inference_summary(use_inds).koff_mean])/60;
        
    kon0_guess = prctile(kon_high_vec,95);
    mcmcOptions.kon0_prior_mean = kon0_guess;
    mcmcOptions.kon0_prior_std = kon0_guess - nanmean(kon_vec); % need a more principled way of setting this
    
    mcmcOptions.kon_hill_prior_mean = -6.4; % from fit to fluo trend
    mcmcOptions.kon_hill_prior_std = 3;
    
    mcmcOptions.koff_hill_prior_mean = -2; % from fit to fluo trend
    mcmcOptions.koff_hill_prior_std = 3;
    
    mcmcOptions.Kd_prior_mean = 5; 
    mcmcOptions.Kd_prior_std = 3;
    
    % assemble into prior vector
    mcmcOptions.prior_mean_vec = repmat(nanmean(mcmcOptions.paramInit,1),size(mcmcOptions.paramInit,1),1);%[mcmcOptions.kon_hill_prior_mean, mcmcOptions.Kd_prior_mean, mcmcOptions.kon0_prior_mean, ...
                                  %mcmcOptions.koff_hill_prior_mean, mcmcOptions.koff0_prior_mean, mcmcOptions.r2_prior_mean];
                                
    mcmcOptions.prior_std_vec = repmat(nanstd(mcmcOptions.paramInit,1),size(mcmcOptions.paramInit,1),1);%[mcmcOptions.kon_hill_prior_std, mcmcOptions.Kd_prior_std, mcmcOptions.kon0_prior_std, ...
                                  %mcmcOptions.koff_hill_prior_std, mcmcOptions.koff0_prior_std, mcmcOptions.r2_prior_std];

%     mcmcOptions.prior_std_vec = repmat(mcmcOptions.prior_std_vec,mcmcOptions.n_chains,1);
    
    % adjust parameter bounds
    paramBounds = mcmcOptions.paramBounds;
    paramBounds(1,1:6) = min(mcmcOptions.prior_mean_vec-2*mcmcOptions.prior_std_vec);
    paramBounds(2,1:6) = max(mcmcOptions.prior_mean_vec+2*mcmcOptions.prior_std_vec);
    
    neg_flags = paramBounds(1,:)<0;
    neg_flags([1 4]) = false;
    paramBounds(1,neg_flags) = 0;
    mcmcOptions.paramBounds = paramBounds;
    