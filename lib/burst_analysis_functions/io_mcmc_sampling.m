function mcmcResults = io_mcmc_sampling(sweepInfo, mcmcOptions)

% make directory to store temporary files

%% extract basic sampling options
n_chains = mcmcOptions.n_chains;
n_mcmc_steps = mcmcOptions.n_mcmc_steps;
prop_factor = mcmcOptions.n_mcmc_prop_size;
param_bounds = mcmcOptions.paramBounds;%(2,1) = -1;
% sweepInfo.n_traces = mcmcOptions.n_traces;

% generate prop sigma vector
prop_sigma_vec = abs(prop_factor*mean(param_bounds));

% initialize parallel pools
sweepInfo.NumWorkers = min([n_chains, sweepInfo.NumWorkers]);
initializePool(sweepInfo)

% initialize stuff for waitbar
% WB = waitbar(0,'conducting mcmc sampling...');
% D = parallel.pool.DataQueue;    
% afterEach(D, @nUpdateWaitbar);
% % 
% % N = sweepInfo.nIterations;
% c = 1;
mcmcResults = struct;

parfor c = 1:n_chains
  
%     sweepInfo = sweepInfo;
  
    % initialize array to store results
    mcmc_param_array = NaN(n_mcmc_steps,length(sweepInfo.paramList));
    
    % initialize array to store logL
    mcmc_loss_array = NaN(n_mcmc_steps,3);

    % store predictions
    ra_prediction_array = NaN(n_mcmc_steps,length(sweepInfo.reactivation_cdf));
    ra_true_prediction_array = NaN(n_mcmc_steps,length(sweepInfo.reactivation_cdf));
    ra_kon_prediction_array = NaN(n_mcmc_steps,length(sweepInfo.time_axis_ra));
    fluo_prediction_array = NaN(n_mcmc_steps,length(sweepInfo.mean_fluo_trend_cm));    
    
    for mcmc_step = 1:n_mcmc_steps
        %%%%%%%%%%%%%%%%%%
        % step 1: generate new step proposal 
        %%%%%%%%%%%%%%%%%%
        % get current guess
        if mcmc_step ~= 1
            current_vals = mcmc_param_array(mcmc_step-1,:);

            % calculate bounds
            lb_array = (param_bounds(1,:)-current_vals)./prop_sigma_vec;
            ub_array = (param_bounds(2,:)-current_vals)./prop_sigma_vec;

            % generate variants
            new_params_prop = current_vals+reshape(prop_sigma_vec'.*trandn(lb_array,ub_array),1,[]);
        else
            new_params_prop = mcmcOptions.paramInit(c,:);%normrnd(mcmcOptions.prior_mean_vec,mcmcOptions.prior_std_vec);
        end
        %%%%%%%%%%%%%%%%%%
        % step 2: run simulation to get likelihood of new params
        %%%%%%%%%%%%%%%%%%
        
        [loss_data, loss_fluo, loss_ra, param_val_vec, pd_fluo_curve, ra_cdf_pd_mean, ra_cdf_pd_mean_true, ra_kon_pd] ...
                                                    = calculate_loss(new_params_prop, sweepInfo, mcmcOptions);
        
        logL_prior = calculate_prior_probabilities(mcmcOptions.prior_mean_vec,mcmcOptions.prior_std_vec, new_params_prop);
        
        loss_total = loss_data + logL_prior;
        
        move_flag = true;
        if mcmc_step > 1
            move_flag = exp(loss_total-mcmc_loss_array(mcmc_step-1,3)) > rand();
        end
        if move_flag 
            mcmc_loss_array(mcmc_step,1) = loss_ra;
            mcmc_loss_array(mcmc_step,2) = loss_fluo;
            mcmc_loss_array(mcmc_step,3) = loss_total;
            
            mcmc_param_array(mcmc_step,:) = param_val_vec;
            fluo_prediction_array(mcmc_step,:) = pd_fluo_curve;
            ra_prediction_array(mcmc_step,:) = ra_cdf_pd_mean;
            ra_true_prediction_array(mcmc_step,:) = ra_cdf_pd_mean_true;
            ra_kon_prediction_array(mcmc_step,:) = ra_kon_pd;
            
        else
            mcmc_loss_array(mcmc_step,:) = mcmc_loss_array(mcmc_step-1,:);            
            mcmc_param_array(mcmc_step,:) = mcmc_param_array(mcmc_step-1,:);
            fluo_prediction_array(mcmc_step,:) = fluo_prediction_array(mcmc_step-1,:);
            ra_prediction_array(mcmc_step,:) = ra_prediction_array(mcmc_step-1,:);
            ra_true_prediction_array(mcmc_step,:) = ra_true_prediction_array(mcmc_step-1,:);
            ra_kon_prediction_array(mcmc_step,:) = ra_kon_prediction_array(mcmc_step-1,:);
            
        end
        % update waitbar
    %     send(D, sweep_step);           
    end
    mcmcResults(c).total_loss = mcmc_loss_array(:,3);
    mcmcResults(c).mcmc_loss_array = mcmc_loss_array;
    mcmcResults(c).mcmc_param_array = mcmc_param_array;
    mcmcResults(c).ra_prediction_array = ra_prediction_array;
    mcmcResults(c).ra_true_prediction_array = ra_true_prediction_array;
    mcmcResults(c).ra_kon_prediction_array = ra_kon_prediction_array;
    mcmcResults(c).fluo_prediction_array = fluo_prediction_array;
end

%%
% delete(WB);

% helper function
% function nUpdateWaitbar(~)
%   waitbar(p/N, WB);
%   p = p + 1;
% end


end