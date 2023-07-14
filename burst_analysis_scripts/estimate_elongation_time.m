clear
close all

addpath(genpath('../lib'));

% attempt to estimate elongation time from autocorrelatio signal
readPath = '..\data\burst_analysis_data\';
load([readPath 'spot_struct.mat'],'spot_struct')

dt = spot_struct(1).tresInterp;
%%
max_time = max([spot_struct.timeInterp]);
timeGrid = 0:dt:max_time;
max_len = length(timeGrid);

% get vector of flags 
use_flag_vec = [spot_struct.useFlag] & [spot_struct.fluoFlag];

% parameters
n_lags = 20;
bootstrap = 1;
trace_weights = ones(size(spot_struct));
n_boots = 100;

% initialize array to store traces
trace_array = zeros(max_len,sum(use_flag_vec));

iter = 1;
for i = find(use_flag_vec)
    ft = ismember(timeGrid,spot_struct(i).timeInterp);
    trace_array(ft,iter) = spot_struct(i).fluoInterp;
    iter = iter + 1;
end    
%%
% run autocorr
[wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors, wt_ddd, ddd_boot_errors] = ...
    weighted_autocorrelation(trace_array, n_lags, bootstrap,n_boots,trace_weights);
  
  