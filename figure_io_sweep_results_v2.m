% Script to generate figures using output of paramers sweeps
clear
close all
addpath(genpath('./lib'))


% point to data directory
dataPath = ['.' filesep 'data' filesep 'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_FINAL';

readPath = [dataPath DateStr suffix filesep];


% point to data directory
%dataPath = ['data' filesep 'burst_analysis_data' filesep];
FolderStr = 'simulation_result';
%readPath = [dataPath FolderStr filesep];

figPath = ['fig' filesep 'burst_analysis_results_reframed' filesep DateStr suffix filesep];
mkdir(figPath)

% load MCMC results
load([readPath 'sweepInfo.mat'],'sweepInfo')
load([readPath 'sweepResults.mat'],'sweepResults')
load([readPath 'mcmcResults.mat'],'mcmcResults')
load([readPath 'mcmcOptions.mat'],'mcmcOptions')

% Concatenate MCMC results
burn_in = 250;
faceAlpha = 0.3;
fluo_norm_factor = 1e4;
bkg_color = [228,221,209]/255;


master_param_array = [];
master_loss_array = [];
fluo_prediction_array = [];
ra_prediction_array = [];
ra_kon_prediction_array = [];
ra_true_prediction_array = [];

for i = 1:length(mcmcResults)
    master_param_array = [master_param_array; mcmcResults(i).mcmc_param_array(burn_in:end,:)];
    master_loss_array = [master_loss_array; mcmcResults(i).mcmc_loss_array(burn_in:end,:)];
    fluo_prediction_array = [fluo_prediction_array; mcmcResults(i).fluo_prediction_array(burn_in:end,:)];
    ra_prediction_array = [ra_prediction_array; mcmcResults(i).ra_prediction_array(burn_in:end,:)];
    
    ra_kon_prediction_array = [ra_kon_prediction_array; mcmcResults(i).ra_kon_prediction_array(burn_in:end,:)];
    ra_true_prediction_array = [ra_true_prediction_array; mcmcResults(i).ra_true_prediction_array(burn_in:end,:)];
end    
% apply normalization factors etc.
master_param_array(:,[3 6]) = master_param_array(:,[3 6])*60; % convert koff and kon to per minute
master_param_array(:,[7 9]) = master_param_array(:,[7 9])/fluo_norm_factor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make figures
close all
n_plot = 25;
light_green = [220 236 203]/256;
dark_green = [122 169 116]/256;
green_array = interp1([1 n_plot],[light_green ; dark_green],1:n_plot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make input Knirps trend figures

rng(235) % for reproducibility
plot_indices = randsample(1:size(sweepInfo.knirps_array_cm,2),n_plot,false);
% WT
wt_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices)
    plot(sweepInfo.time_axis_cm/60,sweepInfo.knirps_array_cm(:,plot_indices(p)),'Color',[green_array(p,:) 0.75],'LineWidth',1)
end    
plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 12])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('[Knirps] (au)')
%grid on
    
%set(gca, 'Color', bkg_color)

wt_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(wt_fig,[figPath 'knirps_wt_trend.png'])
saveas(wt_fig,[figPath 'knirps_wt_trend.pdf'])

% RA

rng(109)
plot_indices_ra = randsample(1:size(sweepInfo.knirps_array_ra,2),n_plot,false);

ra_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices_ra)
    plot(sweepInfo.time_axis_ra/60,sweepInfo.knirps_array_ra(:,plot_indices_ra(p)),'Color',[green_array(p,:) 0.75],'LineWidth',1)
end    
plot(sweepInfo.time_axis_ra(1:2:end)/60,mean(sweepInfo.knirps_array_ra(1:2:end,:),2),'Color','k','LineWidth',1)
scatter(sweepInfo.time_axis_ra(1:2:end)/60,mean(sweepInfo.knirps_array_ra(1:2:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([-10 10])
ylim([0 12])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('[Knirps] (au)')
%grid on
    
%set(gca, 'Color', bkg_color)

ra_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(ra_fig,[figPath 'knirps_ra_trend.png'])
saveas(ra_fig,[figPath 'knirps_ra_trend.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make output trend figs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(243);
err_sig = 1;
lb_prct = 100*normcdf(-err_sig,0,1);
ub_prct = 100*normcdf(err_sig,0,1);
% Plot predicted vs actual IO trends
logL_vec_total = master_loss_array(:,3);
logL_vec_data = master_loss_array(:,1) + master_loss_array(:,2);
% logL_vec_data(master_param_array(:,4)<-1.5) = -Inf;

[logL_vec_data_u, ia, ic] = unique(logL_vec_data);

%% I/O curves 
fluo_mf_vec = mean(fluo_prediction_array,2);
% f999 = prctile(fluo_mf_vec,95);
% f5 = prctile(fluo_mf_vec,5);
% high_ind = find((fluo_mf_vec>=f999)&(master_param_array(:,2)>4.1)&master_param_array(:,3)>3,1);
% 
% low_ind = find((fluo_mf_vec<=f5)&(master_param_array(:,3)<2.5));
% low_ind = low_ind(1);
% [~,si] = sort(logL_vec_data_u,'descend');
n_avg = 25;
opt_params = mean(master_param_array(ia(end-n_avg+1:end),:));
%%%%%%%%%%%%%%%%%
% fluo
% get most likely trend
fluo_pd_opt = mean(fluo_prediction_array(ia(end-n_avg+1:end),:),1)/fluo_norm_factor;
% fluo_pd_high = mean(fluo_prediction_array(high_ind,:),1)/fluo_norm_factor;
% fluo_pd_low = mean(fluo_prediction_array(low_ind,:),1)/fluo_norm_factor;
fluo_pd_mean = mean(fluo_prediction_array)/fluo_norm_factor;
fluo_ub = prctile(fluo_prediction_array,ub_prct,1)/fluo_norm_factor;
fluo_lb = prctile(fluo_prediction_array,lb_prct,1)/fluo_norm_factor;

close all

fluo_io_figure = figure;
cmap = brewermap(8,'Set2');
hold on
fill([sweepInfo.knirps_axis_cm fliplr(sweepInfo.knirps_axis_cm)], [fluo_ub fliplr(fluo_lb)], ...
                                  cmap(2,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)

                         
p1 = errorbar(sweepInfo.knirps_axis_cm, sweepInfo.mean_fluo_trend_cm/fluo_norm_factor,...
                sweepInfo.ste_fluo_trend_cm/fluo_norm_factor,'Color', 'k', 'LineWidth', 1.5,'CapSize',0);
scatter(sweepInfo.knirps_axis_cm, sweepInfo.mean_fluo_trend_cm/fluo_norm_factor,200,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k')                 
p2 = plot(sweepInfo.knirps_axis_cm, fluo_pd_opt, '-.', 'Color', cmap(2,:), 'LineWidth', 3);          

set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
ylabel('mean transcription rate (au)')
%grid on
    
lgd = legend([p1 p2],'experiment','model fit');
xlim([3.2 10.5])
ylim([0 20.55])
%set(gca, 'Color', bkg_color)

fluo_io_figure.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(fluo_io_figure,[figPath 'fluo_io_trend.png'])
saveas(fluo_io_figure,[figPath 'fluo_io_trend.pdf'])
  
delete(lgd)
for n = 1:n_avg
    plot(sweepInfo.knirps_axis_cm, fluo_prediction_array(ia(end-n+1),:)/fluo_norm_factor, '-', 'Color', cmap(2,:), 'LineWidth', 1); 
end

saveas(fluo_io_figure,[figPath 'fluo_io_trend_mean_check.png'])
saveas(fluo_io_figure,[figPath 'fluo_io_trend_mean_check.pdf'])

%%%%%%%%%%%%%%%%%
% reactivation time

% get most likely trend
ra_pd_opt = mean(ra_prediction_array(ia(end-n_avg+1:end),:),1);
% ra_pd_high = mean(ra_prediction_array(high_ind,:),1);
% ra_pd_low = mean(ra_prediction_array(low_ind,:),1);
ra_pd_mean = mean(ra_prediction_array);
ra_ub = prctile(ra_prediction_array,ub_prct,1);
ra_lb = prctile(ra_prediction_array,lb_prct,1);


ra_figure = figure;

hold on
fill([sweepInfo.reactivation_time fliplr(sweepInfo.reactivation_time)]/60, [ra_ub fliplr(ra_lb)], ...
                                  cmap(2,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)

                         
p1 = errorbar(sweepInfo.reactivation_time/60, sweepInfo.reactivation_cdf,...
                sweepInfo.reactivation_cdf_ste,'Color', 'k', 'LineWidth', 1.5,'CapSize',0);
scatter(sweepInfo.reactivation_time/60, sweepInfo.reactivation_cdf,200,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k')              
p2 = plot(sweepInfo.reactivation_time/60, ra_pd_opt, '-.', 'Color', cmap(2,:), 'LineWidth', 3);       

set(gca,'Fontsize',14)

xlabel('time since Knirps export (minutes)')
ylabel('fraction reactivated')
%grid on
    
lgd = legend([p1 p2],'experiment','model fit','Location','southeast');
xlim([0 7])
ylim([0 1])
%set(gca, 'Color', bkg_color)


ra_figure.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(ra_figure,[figPath 'reactivation_trend.png'])
saveas(ra_figure,[figPath 'reactivation_trend.pdf'])

delete(lgd)
for n = 1:n_avg
    plot(sweepInfo.reactivation_time/60, ra_prediction_array(ia(end-n+1),:), '-', 'Color', cmap(2,:), 'LineWidth', 1); 
end
% 
saveas(ra_figure,[figPath 'reactivation_trend_mean_check.png'])
saveas(ra_figure,[figPath 'reactivation_trend_mean_check.pdf'])

%% Make plots of predicted burst parameter IO curves

KD_opt = mean(master_param_array(ia(end-n_avg+1:end),2),1);
KD_opt_koff = mean(master_param_array(ia(end-n_avg+1:end),5),1);
knirps_axis = linspace(0,10);
[~, peg_ind_kon] = min(abs(knirps_axis-KD_opt));
[~, peg_ind_koff] = min(abs(knirps_axis-KD_opt_koff));

% generate array of kon predictions

kon_array = master_param_array(:,3) .* master_param_array(:,2).^-master_param_array(:,1) ./ ...
  (knirps_axis.^-master_param_array(:,1) + master_param_array(:,2).^-master_param_array(:,1));
kon_opt = mean(kon_array(ia(end-n_avg+1:end),:),1);
kon_norm = kon_opt(peg_ind_kon);
kon_pd_mean = mean(kon_array);
kon_ub = prctile(kon_array,ub_prct,1);
kon_lb = prctile(kon_array,lb_prct,1);

koff_array = (master_param_array(:,6) .* master_param_array(:,5).^-master_param_array(:,4) ./ ...
  (knirps_axis.^-master_param_array(:,4) + master_param_array(:,5).^-master_param_array(:,4)));
koff_opt = mean(koff_array(ia(end-n_avg+1:end),:),1);
koff_norm = koff_opt(peg_ind_koff);
koff_pd_mean = mean(koff_array);
koff_ub = prctile(koff_array,ub_prct,1);
koff_lb = prctile(koff_array,lb_prct,1);

r_array = repmat(master_param_array(:,7),1,length(knirps_axis))*3;
r_opt = mean(r_array(ia(end-n_avg+1:end),:),1);
r_norm = r_opt(peg_ind_kon);
r_pd_mean = mean(r_array);
r_ub = prctile(r_array,ub_prct,1);
r_lb = prctile(r_array,lb_prct,1);




cmap_bu = brewermap(8,'Blues');
cmap_gr = brewermap(8,'Greens');
cmap_rd = brewermap(8,'Reds');

close all
% burst_param_figure = figure;
% hold on
% fill([knirps_axis fliplr(knirps_axis)], [kon_ub fliplr(kon_lb)], ...
%                                   cmap_bu(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p1 = plot(knirps_axis, kon_opt, '-.', 'Color', cmap_bu(8,:), 'LineWidth', 3);    
% 
% fill([knirps_axis fliplr(knirps_axis)], [koff_ub fliplr(koff_lb)], ...
%                                   cmap_gr(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p2 = plot(knirps_axis, koff_opt, '-.', 'Color', cmap_gr(8,:), 'LineWidth', 3);  
% 
% ylabel('events per minute')
% ylim([0 6])
% set(gca,'Ytick',[0:1:5])
% 
% yyaxis right 
% 
% fill([knirps_axis fliplr(knirps_axis)], [r_ub fliplr(r_lb)], ...
%                                   cmap_rd(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p3 = plot(knirps_axis, r_opt, '-.', 'Color', cmap_rd(8,:), 'LineWidth', 3);  
% 
% ylim([0 1])
% set(gca,'Ytick',[0:0.2:1])
% ylabel('au per minute')
% set(gca,'Fontsize',14)
% 
% xlabel('[Knirps] (au)')
% grid on
% legend([p1 p2 p3], 'k_{on} (frequency)', 'k_{off} (1/duration)', 'r (amplitude)','Location','northeast')
% set(gca, 'Color', bkg_color)
% 
% burst_param_figure .InvertHardcopy = 'off';
% set(gcf,'color','w');
% ax = gca; 
% ax.YAxis(2).Color = cmap_rd(8,:);
% saveas(burst_param_figure,[figPath 'burst_param_trends.png'])
% saveas(burst_param_figure,[figPath 'burst_param_trends.pdf'])
% 
% % 
% burst_param_figure = figure;
% cmap_bu = brewermap(8,'Blues');
% cmap_gr = brewermap(8,'Greens');
% cmap_rd = brewermap(8,'Reds');
% 
% hold on
% fill([knirps_axis fliplr(knirps_axis)], [kon_ub fliplr(kon_lb)]/kon_norm, ...
%                                   cmap_bu(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p1 = plot(knirps_axis, kon_opt/kon_norm, '-.', 'Color', cmap_bu(8,:), 'LineWidth', 3);    
% 
% fill([knirps_axis fliplr(knirps_axis)], ([koff_ub fliplr(koff_lb)]/koff_norm).^-1, ...
%                                   cmap_gr(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p2 = plot(knirps_axis, (koff_opt/koff_norm).^-1, '-.', 'Color', cmap_gr(8,:), 'LineWidth', 3);  
% 
% 
% fill([knirps_axis fliplr(knirps_axis)], [r_ub fliplr(r_lb)]/r_norm, ...
%                                   cmap_rd(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
% p3 = plot(knirps_axis, r_opt/r_norm, '-.', 'Color', cmap_rd(8,:), 'LineWidth', 3);  
% 
% legend([p1 p2 p3], 'burst frequency (k_{on})', 'burst duration (1/k_{off})', 'burst amplitude (r)','Location','northwest')
% 
% set(gca,'Fontsize',14)
% 
% xlabel('[Knirps] (au)')
% ylabel('normalized parameter magnitude')
% grid on
%     
% xlim([0 8])
% ylim([0 3.5])
% set(gca, 'Color', bkg_color)
% 
% burst_param_figure .InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(burst_param_figure,[figPath 'burst_param_trends_norm.png'])
% saveas(burst_param_figure,[figPath 'burst_param_trends_norm.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kon only plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

% get numeric values for best fit, mean, and CI
kon0_opt = mean(master_param_array(ia(end-n_avg+1:end),3),1)
kon0_mean = mean(master_param_array(:,3),1)
kon0_ste = std(master_param_array(:,3),1)

KD_opt = mean(master_param_array(ia(end-n_avg+1:end),2),1)
KD_mean = mean(master_param_array(:,2),1)
KD_ste = std(master_param_array(:,2),1)

konH_opt = mean(master_param_array(ia(end-n_avg+1:end),1),1)
konH_mean = mean(master_param_array(:,1),1)
konH_ste = std(master_param_array(:,1),1)

% kon_high = kon_array(high_ind,:);
% kon_low = kon_array(low_ind,:);

kon_fig = figure;

hold on
fill([knirps_axis fliplr(knirps_axis)], [kon_ub fliplr(kon_lb)], ...
                                  cmap_bu(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
p1 = plot(knirps_axis, kon_opt, '-.', 'Color', cmap_bu(8,:), 'LineWidth', 3);    


ylim([0 3.5])
xlim([1 9])
ylabel('burst frequency (k_{on})')
set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
%grid on

%set(gca, 'Color', bkg_color)

kon_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% ax = gca; 

pbaspect([3 2 1])

saveas(kon_fig,[figPath 'kon_trend.png'])
saveas(kon_fig,[figPath 'kon_trend.pdf'])

% plot(knirps_axis, kon_high, '-.', 'Color', cmap(8,:), 'LineWidth', 2);
% plot(knirps_axis, kon_low, '-.', 'Color', cmap(7,:), 'LineWidth', 2);
% 
% saveas(kon_fig,[figPath 'kon_trend_alt_trends.png'])
% saveas(kon_fig,[figPath 'kon_trend_alt_trends.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%
% Add binding model fit line
%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 10; % estimated n umber of binding sites from patser screen
n_calc_points = length(kon_opt);

% generate helper functions for fit
% params = [eInteraction, eCoop, kb0, kon0]
kon_pd_fun = @(params) predict_kon_curve_v2(knirps_axis,kon0_opt,params(1),params(2),n_bs,n_calc_points);
loss_fun = @(params) kon_opt-kon_pd_fun(params);

params_fit = lsqnonlin(loss_fun,[1 1], [0 0],[100 100]);
% [kon_curve, n_bound_curve] = predict_kon_curve(repressorCVec,kon0,eInteraction,eCoop,kb0,n_bs,n_points);

kon_binding_pd = kon_pd_fun(params_fit);

p2 = plot(knirps_axis,kon_binding_pd,'Color','k', 'LineWidth', 2);   
legend([p1 p2], 'MCMC inference','binding model fit','Color','w')
uistack(p1,'top')
saveas(kon_fig,[figPath 'kon_trend_binding.png'])
saveas(kon_fig,[figPath 'kon_trend_binding.pdf'])

%%
koff0_opt = mean(master_param_array(ia(end-n_avg+1:end),6),1)
koff0_ste = std(master_param_array(:,6),[],1)

KD_opt = mean(master_param_array(ia(end-n_avg+1:end),5),1)
KD_ste = std(master_param_array(:,5),[],1)

koffH_opt = mean(master_param_array(ia(end-n_avg+1:end),4),1)
konH_ste = std(master_param_array(:,4),[],1)


r_opt = mean(master_param_array(ia(end-n_avg+1:end),7),1)*3
r_ste = std(master_param_array(:,7),[],1)*3

detect_opt = mean(master_param_array(ia(end-n_avg+1:end),9),1)
detect_ste = std(master_param_array(:,9),[],1)

%%
r_fig = figure;

hold on
fill([knirps_axis fliplr(knirps_axis)], [r_ub fliplr(r_lb)], ...
                                  cmap_rd(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
p1 = plot(knirps_axis, r_opt, '-.', 'Color', cmap_rd(8,:), 'LineWidth', 3);    


ylim([0 30])
xlim([1 9])
ylabel('burst amplitude (r)')
set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
grid on

set(gca, 'Color', bkg_color)

r_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(r_fig,[figPath 'r2_trend.png'])
saveas(r_fig,[figPath 'r2_trend.pdf'])



dur_fig = figure;

hold on
fill([knirps_axis fliplr(knirps_axis)], [koff_ub fliplr(koff_lb)].^-1, ...
                                  cmap_gr(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
p1 = plot(knirps_axis, koff_opt.^-1, '-.', 'Color', cmap_gr(8,:), 'LineWidth', 3);    


ylim([0 3])
xlim([1 9])
ylabel('burst duration (1/k_{off})')
set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
grid on

set(gca, 'Color', bkg_color)

dur_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(dur_fig,[figPath 'dur_trend.png'])
saveas(dur_fig,[figPath 'dur_trend.pdf'])

% %% Make bivariate plots (just kon)
% % cmap_gra_short = brewermap(8,'Greys');
% p_dim = 3;
% param_min_vec = [-10 3 2.5];% -5  5  0.75];
% param_max_vec = [ -4 5 4.5];%  0  6.5  0.77];
% 
% % tickCell = {-10:2:-4, 3:5, 3:4};%, -5:2:0, 5:0.5:6.5,[.75 .76 0.77]};
% labelCell = {'H (k_{on})', 'K_D (k_{on})', 'k_{on}^0'};%, 'H (k_{off})', 'k_{off} max','r'};
% labelCell2 = {'H_{on}', 'K_D', 'k_{on}^0'};
% close all
% % colormap_cell = {cmap_bu,cmap_gra_short,cmap_bu,cmap_gr,cmap_gr,cmap_rd};
% 
% grid_size = 75;
% 
% hist_fig = figure;
% cmap_gra = brewermap(128,'Greys');
% cmap_gra(1,:) = bkg_color;
% 
% for i = 1:p_dim
%     param1_i = i;
%     param1_min = floor(10*prctile(master_param_array(:,param1_i),0.1))/10;
%     param1_max = ceil(10*prctile(master_param_array(:,param1_i),99.9))/10;
%     x_axis = linspace(param1_min,param1_max,grid_size);
%     for j = i+1:p_dim
%               
%         param2_i = j;       
%         param2_min = floor(10*prctile(master_param_array(:,param2_i),0.1))/10;
%         param2_max = ceil(10*prctile(master_param_array(:,param2_i),99.9))/10;
%        
%         y_axis = linspace(param2_min,param2_max,grid_size);
%         [x_grid, y_grid] = meshgrid(x_axis,y_axis);
%         [point_density, point_grid,bw] = ksdensity([master_param_array(:,param1_i) master_param_array(:,param2_i)],[x_grid(:),y_grid(:)]);
%         dim = sqrt(length(point_density));
%       
%         % make tile
%         s_ind = (j-1)*p_dim + i;
%         subplot(p_dim,p_dim,s_ind);
%         colormap(cmap_gra);
%         
%         % make contour plot
%         [~,h] = contourf(x_grid, y_grid,reshape(point_density,dim,dim));                                
%         
%         % turn off tick labels for axes not on outer edges
%         if mod(s_ind,p_dim)~=1
%             set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
%         else
% %             set(subplot(p_dim,p_dim,s_ind), 'YTick',tickCell{j})
%             ylabel(labelCell{j})
%         end
%         
%         if s_ind < p_dim*(p_dim-1)+1
%             set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')
%         else
% %             set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
%             xlabel(labelCell{i})
%         end
% %         set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
%     end
%     
%     % make histogram plot
%     s_ind = (i-1)*p_dim + i;
%     subplot(p_dim,p_dim,s_ind);
%     histogram(master_param_array(:,i),x_axis,'FaceColor',cmap_bu(2+i*2,:),'EdgeAlpha',0,'Normalization','probability')
%     hold on
%     % add info about optimal values, uncertainties, etc.
%     param_opt = mean(master_param_array(ia(end-n_avg+1:end),i),1);
%     param_mean = mean(master_param_array(:,i),1);
%     param_ub = prctile(master_param_array(:,i),ub_prct,1);
%     param_lb = prctile(master_param_array(:,i),lb_prct,1);
%                 
%     y_lim = get(gca,'YLim');
%     x_lim = get(gca,'XLim');
%     
%     fill([[param_ub, param_ub] fliplr([param_lb, param_lb])], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], 'k', 'FaceAlpha',0.15,'EdgeAlpha',0)
%     plot([param_mean param_mean], [y_lim(1) y_lim(2)], '-.k', 'LineWidth',1)
%     plot([param_opt param_opt], [y_lim(1) y_lim(2)], '-', 'LineWidth',1.5,'Color',cmap(5,:))
%     
%     % add mean value
%     txt = ['$ \overline{' labelCell2{i} '}$=' num2str(round(param_mean,1)) ];
%     text(0.97*x_lim(1)+0.03*x_lim(2),0.85*y_lim(2),txt,'Interpreter','latex');
%     
%     if i ~= p_dim
%         set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')    
%     else
% %         set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
%         xlabel(labelCell{i})
%     end
%     set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
%     set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
%     
% end    
% hist_fig.InvertHardcopy = 'off';
% hist_fig.Renderer = 'Painters';
% set(gcf,'color','w');
% saveas(hist_fig,[figPath 'mcmc_results_kon.png'])
% saveas(hist_fig,[figPath 'mcmc_results_kon.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make predicted reactivation trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get most likely trend
ra_true_pd_opt = [0 mean(ra_true_prediction_array(ia(end-n_avg+1:end),:),1)];
ra_true_pd_mean = mean(ra_true_prediction_array);
ra_true_ub = [0 prctile(ra_true_prediction_array,ub_prct,1)];
ra_true_lb = [0 prctile(ra_true_prediction_array,lb_prct,1)];

kon_true_pd_opt = mean(ra_kon_prediction_array(ia(end-n_avg+1:end),:),1);
kon_true_pd_mean = mean(ra_kon_prediction_array);
kon_true_ub = prctile(ra_kon_prediction_array,ub_prct,1);
kon_true_lb = prctile(ra_kon_prediction_array,lb_prct,1);


ra_data = [0 sweepInfo.reactivation_cdf];
ra_data_ste = [0 sweepInfo.reactivation_cdf_ste];
ra_axis = [-20 sweepInfo.reactivation_time];
kon_ft = ismember(sweepInfo.time_axis_ra, ra_axis);

ra_figure = figure;

hold on
fill([ra_axis fliplr(ra_axis)]/60, [ra_true_ub fliplr(ra_true_lb)], ...
                                  cmap(5,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)

                         
p1 = errorbar(ra_axis/60, ra_data, ra_data_ste,'Color', 'k', 'LineWidth', 1.5,'CapSize',0);
scatter(ra_axis/60, ra_data,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k')              
p2 = plot(ra_axis/60, ra_true_pd_opt, '-.', 'Color', cmap(5,:), 'LineWidth', 3);       

set(gca,'Fontsize',14)

xlabel('time since Knirps export (minutes)')
ylabel('fraction reactivated')
%grid on
    
xlim([-0.333 5])
ylim([0 1])
set(gca,'Ytick',0:0.25:1)
yyaxis right

fill([ra_axis fliplr(ra_axis)]/60, 60*[kon_true_ub(kon_ft) fliplr(kon_true_lb(kon_ft))], ...
                                  cmap_bu(6,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
                         
p3 = plot(ra_axis/60, kon_true_pd_opt(kon_ft)*60, '--', 'Color', cmap_bu(6,:), 'LineWidth', 3);  
ylabel('burst frequency (events per minute)')

legend([p1 p2 p3],'fraction observed','fraction ON (predicted)','burst frequency (predicted)','Location','southeast')
ax = gca;
ax.YAxis(2).Color = cmap_bu(6,:);
%set(gca, 'Color', bkg_color)
ra_figure.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'Ytick',[0:0.3:1.2])
ylim([0 1.2])

pbaspect([3 2 1])

saveas(ra_figure,[figPath 'reactivation_trend_true.png'])
saveas(ra_figure,[figPath 'reactivation_trend_true.pdf'])

% print HM points to terminal 
ra_axis_hr = linspace(-20,7*20,1e3)/60;

ra_data_hr = interp1(ra_axis/60,ra_data,ra_axis_hr);
[~,mi] = min(abs(ra_data_hr-0.5));
hm_time_data = ra_axis_hr(mi)

ra_pd_hr = interp1(ra_axis/60,ra_true_pd_opt,ra_axis_hr);
[~,mi] = min(abs(ra_pd_hr-0.5));
hm_time_pd = ra_axis_hr(mi)

kon_pd_hr = interp1(ra_axis/60,kon_true_pd_opt(kon_ft)*60,ra_axis_hr);
[~,mi] = min(abs(kon_pd_hr-0.5*max(kon_true_pd_opt(kon_ft)*60)));
hm_time_kon = ra_axis_hr(mi)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make bivariate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap_gra_short = brewermap(8,'Greys');
p_dim = 7;
% param_min_vec = [-10 3 2.5 -5  5  0.75];
% param_max_vec = [ -4 4 4.5  0  6.5  0.77];

% tickCell = {-10:2:-4, 3:5, 3:4, -5:2:0, 5:0.5:6.5,[.75 .76 0.77]};
labelCell = {'H (k_{on})', 'K_D (k_{on})', 'k_{on}^0', 'H (k_{off})', 'K_D (k_{off})','k_{off}^0','r'};
close all
colormap_cell = {cmap_bu,cmap_bu,cmap_bu,cmap_gr,cmap_gr,cmap_gr,cmap_rd};

grid_size = 50;

hist_fig = figure('Position',[100 100 1024 1024]);
cmap_gra = brewermap(128,'Greys');
cmap_gra(1,:) = bkg_color;

for i = 1:p_dim
    param1_i = i;
    param1_min = floor(10*prctile(master_param_array(:,param1_i),0.5))/10;
    param1_max = ceil(10*prctile(master_param_array(:,param1_i),99.5))/10;
    x_axis = linspace(param1_min,param1_max,grid_size);
    for j = i+1:p_dim
              
        param2_i = j;       
        param2_min = floor(10*prctile(master_param_array(:,param2_i),0.5))/10;
        param2_max = ceil(10*prctile(master_param_array(:,param2_i),99.5))/10;
       
        y_axis = linspace(param2_min,param2_max,grid_size);
        [x_grid, y_grid] = meshgrid(x_axis,y_axis);
        [point_density, point_grid,bw] = ksdensity([master_param_array(:,param1_i) master_param_array(:,param2_i)],[x_grid(:),y_grid(:)]);
        dim = sqrt(length(point_density));
      
        % make tile
        s_ind = (j-1)*p_dim + i;
        subplot(p_dim,p_dim,s_ind);
        colormap(cmap_gra);
        [~,h] = contourf(x_grid, y_grid,reshape(point_density,dim,dim));
        
        ax = gca;        
        % turn off tick labels for axes not on outer edges
        if mod(s_ind,p_dim)~=1
            set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
        else
%             set(subplot(p_dim,p_dim,s_ind), 'YTick',tickCell{j})
            ylabel(labelCell{j})
        end
        
        if s_ind < p_dim*(p_dim-1)+1
            set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')
        else
%             set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
            xlabel(labelCell{i})
        end
%         set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
    end
    
    % make histogram plot
    s_ind = (i-1)*p_dim + i;
    subplot(p_dim,p_dim,s_ind);
    histogram(master_param_array(:,i),x_axis,'FaceColor',colormap_cell{i}(6,:),'EdgeAlpha',0,'Normalization','probability')
    
    % add info about optimal values, uncertainties, etc.
    param_opt = mean(master_param_array(ia(end-n_avg+1:end),i),1);
    param_mean = mean(master_param_array(:,i),1);
    param_ub = prctile(master_param_array(:,i),ub_prct,1);
    param_lb = prctile(master_param_array(:,i),lb_prct,1);
    
    hold on
    y_lim = get(gca,'YLim');
    x_lim = get(gca,'XLim');
    
    fill([[param_ub, param_ub] fliplr([param_lb, param_lb])], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], 'k', 'FaceAlpha',0.15,'EdgeAlpha',0)
    plot([param_mean param_mean], [y_lim(1) y_lim(2)], '-.k', 'LineWidth',1)
    plot([param_opt param_opt], [y_lim(1) y_lim(2)], '-', 'LineWidth',1.5,'Color',cmap(5,:))
    
    % add mean value
    txt = ['$' labelCell{i} '$=' num2str(round(param_mean,1)) ];
    text(0.97*x_lim(1)+0.03*x_lim(2),1.1*y_lim(2),txt,'Interpreter','latex');
    
    if i ~= p_dim
        set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')    
    else
%         set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
        xlabel(labelCell{i})
    end
    set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
    set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
    
end    
hist_fig.InvertHardcopy = 'off';
hist_fig.Renderer = 'Painters';
set(gcf,'color','w');
saveas(hist_fig,[figPath 'mcmc_results.png'])
saveas(hist_fig,[figPath 'mcmc_results.pdf'])

%% Look at detection threshold statistics
% close all

kni_axis = linspace(0.5,0.7);

thresh_fig = figure;
cmap_gra_short = brewermap(8,'Greys');
hold on

histogram(master_param_array(:,end),'FaceColor',cmap(1,:),'EdgeAlpha',0,'Normalization','probability')

% add info about optimal values, uncertainties, etc.
param_opt = mean(master_param_array(ia(end-n_avg+1:end),end),1);
param_mean = mean(master_param_array(:,end),1);
param_ub = prctile(master_param_array(:,end),ub_prct,1);
param_lb = prctile(master_param_array(:,end),lb_prct,1);

y_lim = get(gca,'YLim');
x_lim = get(gca,'XLim');

fill([[param_ub, param_ub] fliplr([param_lb, param_lb])], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], 'k', 'FaceAlpha',0.15,'EdgeAlpha',0)
plot([param_mean param_mean], [y_lim(1) y_lim(2)], '-.k', 'LineWidth',1)
plot([param_opt param_opt], [y_lim(1) y_lim(2)], '-', 'LineWidth',1.5,'Color',cmap(5,:))

set(gca,'Fontsize',14)

xlabel('detection threshold (au)')
ylabel('probability')
grid on
set(gca,'Color',bkg_color)    
thresh_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(thresh_fig,[figPath 'threshold_hist.png'])
saveas(thresh_fig,[figPath 'threshold_hist.pdf'])

% generate predicted miss curves
f_axis = sweepInfo.detection_fluo_axis/fluo_norm_factor;
miss_prob_array = master_param_array(:,end).^master_param_array(:,end-1)./...
                (master_param_array(:,end).^master_param_array(:,end-1) + f_axis.^master_param_array(:,end-1));
              
miss_ub = prctile(miss_prob_array,ub_prct,1);
miss_lb = prctile(miss_prob_array,lb_prct,1);
miss_mean = nanmean(miss_prob_array,1);
miss_opt = nanmean(miss_prob_array(ia(end-n_avg+1:end),:),1);

miss_fig = figure;
hold on
fill([f_axis fliplr(f_axis)], [miss_lb fliplr(miss_ub)], ...
                                  cmap_gra_short(3,:), 'FaceAlpha', faceAlpha, 'EdgeColor', 'k','EdgeAlpha',faceAlpha)
plot(f_axis, miss_opt, '-.', 'Color',brighten(cmap_gra_short(3,:),-0.5),'LineWidth',2)

set(gca,'Fontsize',14)

xlabel('MS2 spot intensity (au)')
ylabel('probability of missed detection')
grid on
set(gca,'Color',bkg_color)    
miss_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(miss_fig,[figPath 'miss_prob_plot.png'])
saveas(miss_fig,[figPath 'miss_prob_plot.pdf'])