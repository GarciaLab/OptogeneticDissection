% Script to generate figures using output of paramers sweeps
clear
close all
addpath(genpath('../lib'))

% point to data directory
dataPath = '..\data\burst_analysis_data\';
DateStr = '03-Mar-2022 17-41-41';
suffix = '_v6';
readPath = [dataPath DateStr suffix filesep];

% load data
load([readPath 'sweepResults.mat'],'sweepResults')
load([readPath 'sweepInfo.mat'],'sweepInfo')

ra_weight = 10/11;
cm_weight = 1/11;

lb_full = min(sweepInfo.param_vals);
ub_full = max(sweepInfo.param_vals);
lb = lb_full(:,1:end-2);
% lb(4) = -2;
ub = ub_full(:,1:end-2);
% ub(4) = 0;

% fit_flag = true;%~exist([readPath 'paramFittingResults.mat']);
% if fit_flag
%   
%     fit_struct = struct;
%     n_reps = 25;
%     n_runs = 25;
% %     n_consider = 100;
% 
% %     sweepInfoDS = sweepInfo;
% %     sweepInfoDS.n_traces = 100;
% 
%     myCluster = parcluster('local');
%     max_workers = myCluster.NumWorkers;
%     sweepInfo.NumWorkers = min([n_reps max_workers]);
%     initializePool(sweepInfo);
% 
%     loss_fun = @(params) -calculate_loss(params, sweepInfo);
%     wt_factor = length(sweepInfo.knirps_array_cm) + length(sweepInfo.reactivation_cdf);
%     total_prob_sweep = exp((ra_weight*[sweepResults.ra_r2] + cm_weight*[sweepResults.fluo_r2])*wt_factor);
% %     [maxL, max_i] = max(total_loss_sweep);
% %     [~, si] = sort(total_loss_sweep, 'descend');    
% 
%     % x = fminsearch(loss_fun,params0);
%    
%     % loss_vec = NaN(n_reps,1);
%     % param_array = NaN(n_reps,length(ub));
%     
%     wb = waitbar(0,'Conducting parameter optimization runs');
%     for i = 1:n_reps
%         waitbar(i/n_reps,wb);
%         loss_vec_temp = NaN(n_reps,1);
%         param_array_temp = NaN(n_reps,length(ub));
%         init_indices = randsample(1:length(total_prob_sweep),n_runs,true,total_prob_sweep);
%         init_params = vertcat(sweepResults(init_indices).param_val_vec);
%         for n = 1:n_runs
% 
%             try
%                 [param_fits, total_loss] =...
%                           fmincon(loss_fun, init_params(n,1:end-2),[],[],[],[],lb,ub);
%                 param_array_temp(n,:) = param_fits;
%                 loss_vec_temp(n) = total_loss;        
%             catch
%                 % do nothing
%             end
%         end
%         [loss_vec(i),best_i] = nanmax(loss_vec_temp);
%         param_array(i,:) = param_array_temp(best_i,:);
%         
%         fit_struct(i).loss_vec = loss_vec_temp;
%         fit_struct(i).param_array = param_array_temp;
%         fit_struct(i).init_params = init_params;
%     end
%     delete(wb);
%     save([readPath 'paramFittingResults.mat'], 'param_array')
%     save([readPath 'paramFittingResults.mat'], 'param_array')
%     save([readPath 'paramFittingLoss.mat'], 'loss_vec')
% else
%     disp('Using previous fit results.')
%     load([readPath 'paramFittingResults.mat'])
% %     param_array = paramFittingResults;
%     n_reps = size(param_array,1);
% end    
% Experiment with interpolation
wt_factor = length(sweepInfo.knirps_array_cm) + length(sweepInfo.reactivation_cdf);
total_prob_sweep = exp((ra_weight*[sweepResults.ra_r2] + cm_weight*[sweepResults.fluo_r2])*wt_factor);

pm_array = vertcat(sweepResults.param_val_vec);

% generate array of random values on relevant interval
n_samp = length(total_prob_sweep);
rand_array = rand(n_samp,length(lb_full));
delta_vec = ub_full - lb_full;
pv_array = sweepInfo.param_vals(:,1:5);
% use interpolation to estimate likelihood values
rd_array = rand_array .* delta_vec + lb_full;
% logL_interp = interpn(pm_array,repmat(log(total_prob_sweep'),1,8),rd_array);
logL_interp = interpn(pv_array(:,1),pv_array(:,2),pv_array(:,3),pv_array(:,4),pv_array(:,5),...
                      reshape(log(total_prob_sweep),10,10,10,10,10),...
                      rd_array(:,1),rd_array(:,2),rd_array(:,3),rd_array(:,4),rd_array(:,5),'linear');
d = 5;                    
bw_vec = std(pm_array) .* (4 / ((d+2)*size(pm_array,1))).*(1./(d+4));   
% logL_interp = interpn(param_array_full(:,1),param_array_full(:,2),param_array_full(:,3),param_array_full(:,4),param_array_full(:,5),param_array_full(:,6),param_array_full(:,7),param_array_full(:,8),...
%                       log(total_prob_sweep'),...
%                       rand_param_array(:,1),rand_param_array(:,2),rand_param_array(:,3),rand_param_array(:,4),rand_param_array(:,5),rand_param_array(:,6),rand_param_array(:,7),rand_param_array(:,8));

% resample according to probability
sweep_id_vec = 1:n_samp;
sample_id_vec = randsample(sweep_id_vec,n_samp,true,total_prob_sweep);
param_samp1 = pm_array(sample_id_vec,:);

% use kernel to generate smoothed estimate at likelihood of sample points
% f = mvksdensity(param_samp1(:,1:5),rd_array(:,1:5),'Bandwidth',bw_vec(1:5));

%%
mcmcOptions = struct;
mcmcOptions.n_chains = n_reps; % number of unique chains to run
mcmcOptions.paramBounds = [lb_full ; ub_full] ;
% [~,best_i] = nanmax(samp_weights);
mcmcOptions.paramInit = param_array;
mcmcOptions.n_mcmc_steps = 2.5e3; % number of sampling steps
mcmcOptions.n_mcmc_prop_size = 0.06; % sets scale of proposal jump size
% mcmcOptions.n_init = 10; % top N values from which to select starting point
mcmcOptions.n_traces = 100;
mcmcOptions.burn_in = 100;

cpHMM_path = [dataPath DateStr filesep];
mcmcOptions = getMCMCPriors(mcmcOptions,cpHMM_path);

tic
mcmcResults = io_mcmc_sampling(sweepInfo, mcmcOptions);
toc

save([readPath 'mcmcResults.mat'],'mcmcResults')
save([readPath 'mcmcOptions.mat'],'mcmcOptions')

%%
masterLossArray = vertcat(mcmcResults.mcmc_loss_array);
masterParamArray = vertcat(mcmcResults.mcmc_param_array);
raPredictionArray = vertcat(mcmcResults.ra_prediction_array);
raTruePredictionArray = vertcat(mcmcResults.ra_kon_prediction_array);
fluoPredictionArray = vertcat(mcmcResults.fluo_prediction_array);
jump_frac = mean(diff(masterLossArray(:,3),1,1)~=0)

[maxL, max_i] = nanmax(sum(masterLossArray(:,1:2),2))

figure;
hold on
plot(sweepInfo.mean_fluo_trend_cm)
plot(fluoPredictionArray(max_i,:))
errorbar(mean(fluoPredictionArray,1),std(fluoPredictionArray,1))

figure;
hold on
plot(sweepInfo.reactivation_cdf)
plot(raPredictionArray(max_i,:))
% plot(raTruePredictionArray(max_i,:))
errorbar(mean(raPredictionArray,1),std(raPredictionArray,1))
