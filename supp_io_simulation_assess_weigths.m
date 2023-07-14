% Script to generate figures using output of paramers sweeps
clear
close all
addpath(genpath('../lib'))

% point to data directory
dataPath = ['.' filesep 'data' filesep 'burst_analysis_data' filesep'];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_FINAL';

readPath = [dataPath DateStr suffix filesep];
figPath = ['fig' filesep 'burst_analysis_results_supp' filesep DateStr suffix filesep];
mkdir(figPath)

% load data
load([readPath 'sweepResults.mat'],'sweepResults')
load([readPath 'sweepInfo.mat'],'sweepInfo')

% empirically these weights led to the best visual agreement with
% experimental data trends
n_plot_points = 81;
ra_weight_range = linspace(0.1, 0.9, n_plot_points);%10/11;
fluo_weight_range = 1-ra_weight_range;

% seed random number generator
rng(241);

% get full list of parameter values
param_array_full = vertcat(sweepResults.param_val_vec);
sensitive_koff = param_array_full(:,5)>=3.8;

% initialize array to store optimal parameter values
param_array_opt = NaN(n_plot_points, size(param_array_full,2));
fit_array_opt = NaN(n_plot_points,2);

for n = 1:n_plot_points
    % Use likelihoods to approximate posterior distributio
    total_err_sweep = ra_weight_range(n)*[sweepResults.ra_r2] + fluo_weight_range(n)*[sweepResults.fluo_r2];
    total_err_sweep(sensitive_koff) = -Inf;
    [~,mi] = max(total_err_sweep);
    param_array_opt(n,:) = param_array_full(mi,:);
    fit_array_opt(n,1) = sweepResults(mi).ra_r2;
    fit_array_opt(n,2) = sweepResults(mi).fluo_r2;
end

% Get I/O predictions for top performers
pd_fluo_array = NaN(n_plot_points, length(sweepInfo.knirps_axis_cm));
pd_ra_array = NaN(n_plot_points, length(sweepInfo.reactivation_time));

tic
for n = 1:n_plot_points
    params_prop = param_array_opt(n,:);

    % conduct WT simulations
    [params_prop, ~, ~, pd_fluo_array(n,:), ~] = io_prediction_cm(params_prop, sweepInfo);
    
    % conduct RA simulations
    [~, ~, pd_ra_array(n,:), ra_cdf_pd_true, kon_curve_ra] = io_prediction_ra(params_prop, sweepInfo);
end
toc

%% 
rng(27);
plot_ind = 66;
n_reps = 5;
point_fluo_array = NaN(n_reps, length(sweepInfo.knirps_axis_cm));
point_ra_array = NaN(n_reps, length(sweepInfo.reactivation_time));

wt_factor = length(sweepInfo.knirps_array_cm) + length(sweepInfo.reactivation_cdf);
total_prob_sweep = exp((ra_weight_range(plot_ind)*[sweepResults.ra_r2] + fluo_weight_range(plot_ind)*[sweepResults.fluo_r2])*wt_factor);
sample_inds = randsample(1:length(total_prob_sweep),n_reps,true,total_prob_sweep);

tic
for n = 1:n_reps
    params_prop = param_array_full(sample_inds(n),:);
    % conduct WT simulations
    [params_prop, ~, ~, point_fluo_array(n,:), ~] = io_prediction_cm(params_prop, sweepInfo);
    
    % conduct RA simulations
    [~, ~, point_ra_array(n,:), ra_cdf_pd_true, kon_curve_ra] = io_prediction_ra(params_prop, sweepInfo);
end
toc

% Generate figures
plot_range = 1:n_plot_points;

f_norm_factor = 1e4;
close all

fluo_fig = figure;
cmap = brewermap(n_plot_points,'Spectral');
hold on
colormap(cmap)

iter = 1;
for n = plot_range
    plot(sweepInfo.knirps_axis_cm,pd_fluo_array(n,:)/f_norm_factor,'-','Color',[cmap(n,:) 0.5],'LineWidth',2);
    iter = iter + 1;
end
cb = colorbar;
xlim([3 10])
ylim([0 25])
p2 = plot(-1,-1,'-.','Color','k','LineWidth',2);
p1 = plot(sweepInfo.knirps_axis_cm,sweepInfo.mean_fluo_trend_cm/f_norm_factor,'Color','k','LineWidth',5);
p3 = plot(sweepInfo.knirps_axis_cm,nanmean(point_fluo_array,1)/f_norm_factor,'-','Color',[1 1 1]*0.66,'LineWidth',5);

legend([p1 p2 p3],'experimental data','simulation predictions','selected weights')
clim([ra_weight_range(1) ra_weight_range(end)])
set(gca,'FontSize',14)
ylabel(cb,'reactivation weight (wra)')
xlabel('[Knirps] (au)')
ylabel('mean transcription rate (au)')
saveas(fluo_fig, [figPath 'fluo_io_range.png'])
saveas(fluo_fig, [figPath 'fluo_io_range.pdf'])

% set(cb,'xticklabels',plot_range)
% RA
ra_fig = figure;
colormap(cmap)
hold on
iter = 1;
for n = plot_range
    plot(sweepInfo.reactivation_time/60,pd_ra_array(n,:),'-','Color',[cmap(n,:) 0.5],'LineWidth',1.5);
    iter = iter + 1;
end
p1 = plot(sweepInfo.reactivation_time/60,sweepInfo.reactivation_cdf,'Color','k','LIneWidth',5);
p3 = plot(sweepInfo.reactivation_time/60,mean(point_ra_array,1),'-','Color',[1 1 1]*0.66,'LineWidth',5);
cb = colorbar;
clim([ra_weight_range(1) ra_weight_range(end)])
xlim([0 7])
ylim([0 1.05])
p2 = plot(-1,-1,'-.','Color',[1 1 1]/3,'LineWidth',2);

set(gca,'FontSize',14)
ylabel(cb,'reactivation weight (wra)')
xlabel('time (min) relative to perturbation (minutes)')
ylabel('fraction reactivated')
legend([p1 p2 p3],'experimental data','simulation predictions','selected weights','Location','southeast')

saveas(fluo_fig, [figPath 'ra_range.png'])
saveas(fluo_fig, [figPath 'ra_range.pdf'])