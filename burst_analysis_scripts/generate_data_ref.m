clear 
close all

% specify key hyperparameters
data_ref = struct;
data_ref.knirps_norm_factor = 1e5;
data_ref.knirps_offset = 3.75e5 / data_ref.knirps_norm_factor;
data_ref.min_dp = 10;
data_ref.ap_bounds = [-0.02 0.02];
data_ref.knirps_bins_cpHMM = [0 3.75 5.25 10];
data_ref.fluo_io_min_time = 5; 
data_ref.knirps_cal_slope = 1.243;
data_ref.knirps_cal_intercept = 1.079e5 / data_ref.knirps_norm_factor; %NL: dividng everything through by 1e5 for simplicity
save('data_ref')