% Script to generate figures using output of paramers sweeps
clear
close all
addpath(genpath('./lib'))

% point to data directory
dataPath = ['data' filesep 'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';
suffix_no = '_no_koff';
suffix_yes = '_v3';
readPath_no = [dataPath DateStr suffix_no filesep];
readPath_yes = [dataPath DateStr suffix_yes filesep];

% generate figure path
figPath = ['.\fig\burst_analysis_results\supp_figs\'];
mkdir(figPath)

% load datasets
load([readPath_no 'sweepResults.mat'],'sweepResults')
sweepResultsNo = sweepResults;
clear sweepResults

load([readPath_yes 'sweepResults.mat'],'sweepResults')
sweepResultsYes = sweepResults;
clear sweepResults

% extract loss vectors
cm_loss_vec_no_koff = [sweepResultsNo.fluo_io_fit];
ra_loss_vec_no_koff = [sweepResultsNo.ra_fit];
total_loss_vec_no_koff = cm_loss_vec_no_koff + ra_loss_vec_no_koff;
[~,si_no] = sort(total_loss_vec_no_koff,'descend');

cm_loss_vec = [sweepResultsYes.fluo_io_fit];
ra_loss_vec = [sweepResultsYes.ra_fit];
total_loss_vec = cm_loss_vec+ ra_loss_vec;
[~,si_yes] = sort(total_loss_vec,'descend');
%% Make figures
n_plot = 5e3;
bkg_color = [228,221,209]/255;

scatter_fig = figure;
cmap = brewermap([],'Set2');
hold on

scatter(ra_loss_vec(si_yes(1:n_plot)),cm_loss_vec(si_yes(1:n_plot)),20,'MarkerFaceColor',cmap(2,:),'MarkerFaceAlpha',0.4,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.4)
scatter(ra_loss_vec_no_koff(si_no(1:n_plot)),cm_loss_vec_no_koff(si_no(1:n_plot)),20,'MarkerFaceColor','k','MarkerFaceAlpha',1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2)
xlabel('reactivation fit error')
ylabel('mean activity fit error')
grid on
set(gca, 'Color', bkg_color)
xlim([-2 0])
ylim([-2 0])

saveas(scatter_fig,[figPath 'loss_scatter_fit.png'])
saveas(scatter_fig,[figPath 'loss_scatter_fit.pdf'])
