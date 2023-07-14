function sweepInfo = getMarkovSystemInfo_v2(sweepInfo,cpHMMPath)

    % load inference results
    load([cpHMMPath filesep 'inference_summary.mat'],'inference_summary')

    sweepInfo.memory = 7;
    sweepInfo.deltaT = 20;
    sweepInfo.t_MS2 = 1.4;        
    
    knirps_vec = [inference_summary.knirps_mean];
    knirps_ste_vec = [inference_summary.knirps_ste];
    knirps_axis = linspace(0,1.1*max(knirps_vec));
    nBoots = 100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OFF RATE%
    % when we're not sweeping koff slope, use mean and quartiles as bounds
    koff_vec = [inference_summary.koff_mean]/60;
    koff_ste_vec = [inference_summary.koff_ste]/60;
    
    if ~sweepInfo.koff_slope_sweep_flag        
        koff0_lower = prctile(koff_vec,25);
        koff0_upper = prctile(koff_vec,75);
        koff0_guess = nanmean(koff_vec);
        koff0_sigma = nanstd(koff_vec);
    else
        % otherwise, conduct simple bootstrapped linear fit to obtain
        % bounds on intercept and slope
      
        nan_ft = ~isnan(koff_vec);
        koff_weight_vec = (koff_ste_vec + 0.01*nanmean(koff_vec)).^-2; % extra factor prevents over-fitting of low-variance points        
        [koff_trend, ~, ~, ~, param_95, param_05, koff_array]...
                                              = fit_lin_trend(...
                                              knirps_axis,knirps_vec(nan_ft),knirps_ste_vec(nan_ft), koff_vec(nan_ft),...
                                                koff_ste_vec(nan_ft),koff_weight_vec(nan_ft).^2,nBoots);        
        koff0_upper = param_95(1);
        koff0_lower = param_05(1);
        koff0_sigma = std(koff_array,[],2);
        koff0_sigma = koff0_sigma(1);
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
    [kon_trend, ~, ~, ~, kon_95, kon_05, kon_array]...
                                          = fit_lin_trend(...
                                          knirps_axis,knirps_vec(nan_ft),knirps_ste_vec(nan_ft), kon_vec(nan_ft),...
                                            kon_ste_vec(nan_ft),kon_weight_vec(nan_ft).^2,nBoots);        
    kon0_upper = kon_95(1);
    kon0_lower = kon_05(1);
    kon0_sigma = std(kon_array,[],2);
    kon0_sigma = kon0_sigma(1);
    kon0_guess = kon_trend(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initiation rate bounds
    r1_mean = nanmean([inference_summary.r1_mean]);
    r2_guess = nanmean([inference_summary.r2_mean]);
    r2_sigma = nanstd([inference_summary.r2_mean]);
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
    
    % generate array of defaultvalues
    % [HC_kon KD_kon kon0 HC_koff KD_koff koff0 r2 HC_detect KD_detect]
    sweepInfo.defaultValues = [-5 3.5 kon0_guess -2 3.5 koff0_guess r2_guess*sweepInfo.deltaT NaN NaN];
    sweepInfo.priorSigmas = [3 1 kon0_sigma 2 3.5 koff0_sigma r2_sigma*sweepInfo.deltaT NaN NaN];