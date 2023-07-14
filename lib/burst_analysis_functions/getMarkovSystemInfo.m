function sweepInfo = getMarkovSystemInfo(sweepInfo,cpHMMPath)

    % load inference results
    load([cpHMMPath filesep 'inference_summary.mat'],'inference_summary')

    sweepInfo.memory = 7;
    sweepInfo.deltaT = 20;
    sweepInfo.t_MS2 = 1.4;
    
    % get average parameter values
    exp_type_cell = {inference_summary.exp_type};
    use_inds = contains(exp_type_cell,{'Wild','ON'});
    high_ind = contains(exp_type_cell,{'HIGH'});
    
    % calculate bounding parameters for sweep    
    kon0_lower = prctile([inference_summary(use_inds).kon_mean],75)/60; % used as lower bound for kon0
    kon_high_vec = inference_summary(high_ind).kon_mean/60; % used as "best guess"
    kon0_upper = nanmax([inference_summary(use_inds).kon_mean])/60 * 1.1; % used to set upper bound
    
    % when we're not sweeping koff slope, use mean and quartiles as bounds
    koff_vec = [inference_summary(use_inds).koff_mean]/60;
    koff_ste_vec = [inference_summary(use_inds).koff_ste]/60;
    koff0_guess = nanmean(koff_vec);
    if ~sweepInfo.koff_slope_sweep_flag        
        koff0_lower = prctile(koff_vec,25);
        koff0_upper = prctile(koff_vec,75);
        koff_slope_upper = 0;
        koff_slope_lower = 0;
    else
        % otherwise, conduct simple bootstrapped linear fit to obtain
        % bounds on intercept and slope
        knirps_vec = [inference_summary(use_inds).knirps_mean];
        knirps_ste_vec = [inference_summary(use_inds).knirps_ste];
        knirps_axis = linspace(0,1.1*max(knirps_vec));
        nan_ft = ~isnan(koff_vec);
        koff_weight_vec = (koff_ste_vec + 0.01*nanmean(koff_vec)).^-2;
        nBoots = 100;
        [~, ~, ~, ~, param_95, param_05, ~]...
                                              = fit_lin_trend(...
                                              knirps_axis,knirps_vec(nan_ft),knirps_ste_vec(nan_ft), koff_vec(nan_ft),...
                                                koff_ste_vec(nan_ft),koff_weight_vec(nan_ft).^2,nBoots);
        
        koff0_upper = param_95(1);
        koff0_lower = param_05(1);
        
        koff_slope_upper = param_95(2);
        koff_slope_lower = param_05(2);
    end
    % initiation rate bounds
    r1_mean = nanmean([inference_summary(use_inds).r1_mean]);
    r2_mean = nanmean([inference_summary(use_inds).r2_mean]);
    r2_lower = prctile([inference_summary(use_inds).r2_mean],25);
    r2_upper = prctile([inference_summary(use_inds).r2_mean],75);
    
    % noise in MS2 signal
    noise_mean = nanmedian([inference_summary(use_inds).noise_mean]);
    
    % calculate/stipulate prior distributions over parameter values
    sweepInfo.r2_prior_mean = nanmean([inference_summary(use_inds).r2_mean])*sweepInfo.deltaT;
    sweepInfo.r2_prior_std = nanstd([inference_summary(use_inds).r2_mean])*sweepInfo.deltaT;

    sweepInfo.koff_prior_mean = nanmean([inference_summary(use_inds).koff_mean])/60;
    sweepInfo.koff_prior_std = nanstd([inference_summary(use_inds).koff_mean])/60;
        
    kon0_guess = prctile(kon_high_vec,95);
    sweepInfo.kon_prior_mean = kon0_guess;
    sweepInfo.kon_prior_std = nanstd(kon_high_vec); % need a more principled way of setting this
    
    sweepInfo.kon_hill_prior_mean = 6.4; % from fit to fluo trend
    sweepInfo.kon_hill_prior_std = 1.5;
    
    sweepInfo.Kd_prior_mean = 5; 
    sweepInfo.Kd_prior_std = 2;
    
    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    sweepInfo.R2_orig = [-kon0_guess  koff0_guess; 
                        kon0_guess -koff0_guess];    
    
    sweepInfo.r_emission = [r1_mean r2_mean]*sweepInfo.deltaT;
    sweepInfo.pi0 = [koff0_guess/(kon0_guess +koff0_guess) kon0_guess /(kon0_guess +  koff0_guess)];    
    sweepInfo.noise = noise_mean;
    
    sweepInfo.kon0_lower = kon0_lower;
    sweepInfo.kon0_upper = kon0_upper;
    
    sweepInfo.koff0_lower = koff0_lower;
    sweepInfo.koff0_upper = koff0_upper;
    sweepInfo.koff_mean = koff0_guess;
    sweepInfo.koff_slope_upper = koff_slope_upper;
    sweepInfo.koff_slope_lower = koff_slope_lower;
    
    sweepInfo.r2_lower = r2_lower*sweepInfo.deltaT;
    sweepInfo.r2_mean = r2_mean*sweepInfo.deltaT;
    sweepInfo.r2_upper = r2_upper*sweepInfo.deltaT;
    