function [loss_curr, fluo_loss, ra_loss, params_prop, pd_fluo_curve, ra_cdf_pd_mean, ...
                         ra_cdf_pd_true, kon_curve_ra] = calculate_loss(params_prop, sweepInfo, mcmcOptions)

    % add dummy slots for thresholding parameters 
    paramList = sweepInfo.paramList;
    if length(params_prop) < length(paramList)
        detect_flags = contains(paramList,'detect');    
        param_val_vec_temp = NaN(size(paramList));
        param_val_vec_temp(~detect_flags) = params_prop;    
        params_prop = param_val_vec_temp;
    end
    
    % conduct WT simulations
    [params_prop, sweepInfo, ~, pd_fluo_curve, ~] = io_prediction_cm(params_prop, sweepInfo);

    % conduct RA simulations
    [~, ~, ra_cdf_pd_mean, ra_cdf_pd_true, kon_curve_ra] = io_prediction_ra(params_prop, sweepInfo);       

    %%%%%%%%%%%%%%%%%%%%
    % step 3: Do MH move
    %%%%%%%%%%%%%%%%%%%%
    % likelihood given data
    % note that normalization by mean is equivalent to assuming a Gaussian
    % model with sigma = <mean of empirical data>/2
    fluo_diffs = (pd_fluo_curve'-sweepInfo.mean_fluo_trend_cm) / nanmean(sweepInfo.mean_fluo_trend_cm);
    cdf_diffs = (ra_cdf_pd_mean-sweepInfo.reactivation_cdf) / nanmean(sweepInfo.reactivation_cdf);
    % NL: I've observed that upweighting RA leads to optimal joint fitting
    fluo_loss = -nanmean(fluo_diffs.^2)*mcmcOptions.fluo_weight;
    ra_loss = -nanmean(cdf_diffs.^2)*mcmcOptions.ra_weight;
    loss_curr =  (fluo_loss + ra_loss)*(length(cdf_diffs)+length(fluo_diffs));