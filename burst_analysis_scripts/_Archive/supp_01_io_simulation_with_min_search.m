% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath(['..' filesep 'lib']))

% point to data directory
dataPath = ['..' filesep 'data' filesep 'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_ignore';
outPath = [dataPath DateStr suffix filesep];
mkdir(outPath)

% inisitalize sweep structure and set some parameters
sweepInfo.nParamIncrement = 5;
sweepInfo.keep_prediction_flag = false;
sweepInfo.r2_sweep_flag = false;
sweepInfo.koff_sweep_flag = true;
sweepInfo.koff_slope_sweep_flag = true;

% set list of parameters to sample   
sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_koff','koff0','r2','HC_detect','KD_detect'};
sweepInfo.simFlags = [true,   true,    true,    false(1,5)];
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'koff0')) = sweepInfo.koff_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'HC_koff')) = sweepInfo.koff_slope_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'r2')) = sweepInfo.r2_sweep_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sweep info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call wrapper script;
sweepInfo = initialize_io_sweep(dataPath, sweepInfo, DateStr);                                    

% sweepInfo.nFit = length(sweepInfo.paramList);
sweepInfo.nSim = sum(sweepInfo.simFlags);
sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nSim;

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);
sweepResults = initializeSweepValues(sweepInfo, sweepResults, []);              
    
tic
sweepResults = sweep_loop_v3(sweepInfo,sweepResults);    
toc

%% attempt to use built-in function to search for minimum
sweepInfoDS = sweepInfo;
sweepInfoDS.n_traces = 100;
loss_fun = @(params) calculate_loss(params, sweepInfo);
total_loss_sweep = [sweepResults.ra_r2] + [sweepResults.fluo_r2];
[maxL, max_i] = max(total_loss_sweep);
[~, si] = sort(total_loss_sweep, 'descend');
best_params = sweepResults(max_i).param_val_vec;

% x = fminsearch(loss_fun,params0);
lb_full = min(sweepInfo.param_vals);
ub_full = max(sweepInfo.param_vals);
ub_full([1 4]) = -1;
lb_full(4) = -3;
lb = lb_full(:,1:end-2);
ub = ub_full(:,1:end-2);

n_runs = 50;
loss_vec = NaN(n_runs,1);
param_array = NaN(n_runs,length(ub));
fit_struct = struct;
parfor n = 1:n_runs
%     samp_index = 1:length(sweepResults);
%     samp_i = randsample(samp_index,1);
    try
        [param_fits, total_loss,fit_struct(n).exitflag,fit_struct(n).output,...
            fit_struct(n).lambda,fit_struct(n).grad,fit_struct(n).hessian] =...
                  fmincon(loss_fun, best_params(1:end-2),[],[],[],[],lb,ub);
        param_array(n,:) = param_fits;
        loss_vec(n) = total_loss;        
    catch
        % do nothing
    end
end
%% Use MCMC to estimate parameter uncertainties
samp_weights = abs(1./loss_vec);

mcmcOptions = struct;
mcmcOptions.n_chains = 24; % number of unique chains to run
mcmcOptions.paramBounds = [lb_full ; ub_full] ;
[~,best_i] = nanmax(samp_weights);
mcmcOptions.paramInit = repmat(param_array(best_i,:),mcmcOptions.n_chains,1);
mcmcOptions.n_mcmc_steps = 5e2; % number of samplinqg steps
mcmcOptions.n_mcmc_prop_size = 0.04; % sets scale of proposal jump size
% mcmcOptions.n_init = 10; % top N values from which to select starting point
mcmcOptions.n_traces = 50;
mcmcOptions.burn_in = 100;

cpHMM_path = [dataPath DateStr filesep];
mcmcOptions = getMCMCPriors(mcmcOptions,cpHMM_path);

tic
mcmcResults = io_mcmc_sampling(sweepInfoDS, mcmcOptions);
toc

save([outPath 'mcmcResults.mat'],'mcmcResults')
save([outPath 'mcmcOptions.mat'],'mcmcOptions')

%%
test = vertcat(mcmcResults.mcmc_param_array);