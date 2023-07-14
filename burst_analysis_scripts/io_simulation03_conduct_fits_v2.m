% Script to generate figures using output of paramers sweeps
clear
close all
addpath(genpath('../lib'))

% point to data directory
dataPath = ['..' filesep 'data' filesep 'burst_analysis_data' filesep'];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_FINAL';
% suffix = '_RR2';
readPath = [dataPath DateStr suffix filesep];

% load data
load([readPath 'sweepResults.mat'],'sweepResults')
load([readPath 'sweepInfo.mat'],'sweepInfo')

% empirically these weights led to the best visual agreement with
% experimental data trends
ra_weight = 3/4;
fluo_weight = 1/4;
% ra_weight = 1/2;
% fluo_weight = 1/2;

lb_full = min(sweepInfo.param_vals);
ub_full = max(sweepInfo.param_vals);
lb = lb_full(:,1:end-2);
% lb(4) = -2;
ub = ub_full(:,1:end-2);
% ub(4) = 0;

% seed random number generator
rng(247);

% Use likelihoods to approximate posterior distribution
% resample to estimate priors for key sweep parameters
wt_factor = length(sweepInfo.knirps_array_cm) + length(sweepInfo.reactivation_cdf);
total_prob_sweep = exp((ra_weight*[sweepResults.ra_r2] + fluo_weight*[sweepResults.fluo_r2])*wt_factor);

%% get full list of parameter values
param_array_full = vertcat(sweepResults.param_val_vec);

% generate array of random values on relevant interval
n_samp = length(total_prob_sweep);

% resample according to probability
sweep_id_vec = 1:n_samp;
sample_id_vec = randsample(sweep_id_vec,n_samp,true,total_prob_sweep);
param_array_resamp = param_array_full(sample_id_vec,:);

n_chains = 24;
mcmcOptions = struct;
mcmcOptions.ra_weight = ra_weight;
mcmcOptions.fluo_weight = fluo_weight;
mcmcOptions.n_chains = n_chains; % number of unique chains to run
mcmcOptions.paramBounds = [lb_full ; ub_full] ;
if isfield(sweepInfo,'kon_only_flag')
    mcmcOptions.kon_only_flag = sweepInfo.kon_only_flag;
else    
    mcmcOptions.kon_only_flag = false;
end    
% [~,best_i] = nanmax(samp_weights);
% mcmcOptions.paramInit = param_array;
mcmcOptions.n_mcmc_steps = 2.5e3; % number of sampling steps
mcmcOptions.n_mcmc_prop_size = 0.05; % sets scale of proposal jump size
% mcmcOptions.n_init = 10; % top N values from which to select starting point
sweepInfo = rmfield(sweepInfo,'n_traces');
sweepInfo.n_traces_ra = 215;
sweepInfo.n_traces_cm = 621;
mcmcOptions.burn_in = 100;

% cpHMM_path = [dataPath DateStr filesep];
mcmcOptions = getMCMCPriors_v2(mcmcOptions,sweepInfo,param_array_resamp);

tic
mcmcResults = io_mcmc_sampling(sweepInfo, mcmcOptions);
toc

save([readPath 'mcmcResults.mat'],'mcmcResults')
save([readPath 'mcmcOptions.mat'],'mcmcOptions')

%%
% masterLossArray = vertcat(mcmcResults.mcmc_loss_array);
% masterParamArray = vertcat(mcmcResults.mcmc_param_array);
% raPredictionArray = vertcat(mcmcResults.ra_prediction_array);
% raTruePredictionArray = vertcat(mcmcResults.ra_kon_prediction_array);
% fluoPredictionArray = vertcat(mcmcResults.fluo_prediction_array);
% jump_frac = mean(diff(masterLossArray(:,3),1,1)~=0)
% 
% [maxL, max_i] = nanmax(sum(masterLossArray(:,1:2),2))
% 
% figure;
% hold on
% plot(sweepInfo.mean_fluo_trend_cm)
% plot(fluoPredictionArray(max_i,:))
% errorbar(mean(fluoPredictionArray,1),std(fluoPredictionArray,1))
% 
% figure;
% hold on
% plot(sweepInfo.reactivation_cdf)
% plot(raPredictionArray(max_i,:))
% % plot(raTruePredictionArray(max_i,:))
% errorbar(mean(raPredictionArray,1),std(raPredictionArray,1))
