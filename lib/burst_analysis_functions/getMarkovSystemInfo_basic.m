function sweepInfo = getMarkovSystemInfo_basic(sweepInfo,cpHMMPath)

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
    koff_mean = nanmean(koff_vec);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ON RATE%
    % when we're not sweeping koff slope, use mean and quartiles as bounds
    kon_vec = [inference_summary.kon_mean]/60;
    kon_mean = nanmean(kon_vec);
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initiation rate bounds
    r1_mean = nanmean([inference_summary.r1_mean]);
    r2_mean = nanmean([inference_summary.r2_mean]);  
    
    % noise in MS2 signal
    noise_mean = nanmedian([inference_summary.noise_mean]);        
    
    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    sweepInfo.R2_orig = [-kon_mean  koff_mean; 
                          kon_mean -koff_mean];    
    
    sweepInfo.r_emission = [r1_mean r2_mean]*sweepInfo.deltaT;
    sweepInfo.pi0 = [koff_mean/(kon_mean +koff_mean) kon_mean /(kon_mean +  koff_mean)];    
    sweepInfo.noise = noise_mean;
    