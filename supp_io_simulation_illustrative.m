% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../lib'))

% point to data directory
dataPath = ['.' filesep 'data' filesep 'burst_analysis_data' filesep];
DateStr = '03-Mar-2022 17-41-41';
suffix = '_FINAL_supp';
outPath = [dataPath DateStr suffix filesep];
mkdir(outPath)

figPath = ['fig' filesep 'burst_analysis_results_supp' filesep DateStr suffix filesep];
mkdir(figPath)

% seed random number generator
rng(237);

% inisitalize sweep structure and set some parameters
sweepInfo.nParamIncrement = 1;
sweepInfo.keep_prediction_flag = true;
sweepInfo.r2_sweep_flag = false;
sweepInfo.kon_sweep_flag = false; 
sweepInfo.koff_sweep_flag = false;
if contains(suffix,'kon_only')
    sweepInfo.koff_kd_sweep_flag = false;
    sweepInfo.koff_slope_sweep_flag = false;
else
    sweepInfo.koff_kd_sweep_flag = true;
    sweepInfo.koff_slope_sweep_flag = true;
end
% set list of parameters to sample   
sweepInfo.paramList = {'HC_kon','KD_kon','kon0','HC_koff','KD_koff','koff0','r2','HC_detect','KD_detect'};
sweepInfo.simFlags = [true,   true,   false(1,7)];
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'koff0')) = sweepInfo.koff_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'HC_koff')) = sweepInfo.koff_slope_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'KD_koff')) = sweepInfo.koff_kd_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'kon0')) = sweepInfo.kon_sweep_flag;
sweepInfo.simFlags(strcmp(sweepInfo.paramList,'r2')) = sweepInfo.r2_sweep_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sweep info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call wrapper script;
sweepInfo = initialize_io_sweep(dataPath, sweepInfo, DateStr);                                    

% sweepInfo.nFit = length(sweepInfo.paramList);
sweepInfo.nSim = sum(sweepInfo.simFlags);
sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nSim;
sweepInfo.kon_only_flag = false;
if ~sweepInfo.koff_slope_sweep_flag && ~sweepInfo.koff_kd_sweep_flag
    sweepInfo.kon_only_flag = true;
end

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);

rng(234);
param_val_vec = [-6.05 3.8 (2.80/60) -3.19 3.33 (5.85/60) 7.18e4 NaN NaN]; % optimal values returned by MCMC sampling
sweepResults = initializeSweepValues(sweepInfo, sweepResults, param_val_vec);              
   
tic
sweepResults = sweep_loop_v3(sweepInfo,sweepResults);    
toc
%% Make figures
close all
n_plot = 3;
light_green = [220 236 203]/256;
dark_green = [122 169 116]/256;
green_array = interp1([1 n_plot],[light_green ; dark_green],1:n_plot);

plot_indices = 12;
burst_param_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices)
    plot(sweepInfo.time_axis_cm/60,sweepResults(1).knirps_traces_cm(:,plot_indices(p)),'Color',[dark_green 1],'LineWidth',3)
end    
% plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
% scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 12])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('[Knirps] (au)')
%grid on
    
%set(gca, 'Color', bkg_color)

burst_param_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(burst_param_fig,[figPath 'knirps_wt_trend.png'])
saveas(burst_param_fig,[figPath 'knirps_wt_trend.pdf'])

% burst parameter trends

burst_param_fig = figure('Position',[100 100 512 256]);
cmap = brewermap([],'Set2');
hold on
for p = 1:length(plot_indices)
    plot(sweepInfo.time_axis_cm/60,sweepResults(1).kon_curves_cm(:,plot_indices(p))*60,'Color',cmap(3,:),'LineWidth',3)
    plot(sweepInfo.time_axis_cm/60,(sweepResults(1).koff_curves_cm(:,plot_indices(p))*60),'Color',cmap(2,:),'LineWidth',3)
end    
% plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
% scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 6])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('transitions per minute')
%grid on
    
%set(gca, 'Color', bkg_color)

burst_param_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(burst_param_fig,[figPath 'parameter_wt_trend.png'])
saveas(burst_param_fig,[figPath 'parameter_wt_trend.pdf'])

% PROMOTER
pm_color = [115 142 192]/256;
promoter_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices)
    stairs(sweepResults(1).time_fine_cm/60,sweepResults(1).promoter_state_array_cm(:,plot_indices(p))-1,'Color',pm_color,'LineWidth',2)
end    
% plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
% scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 1.1])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('promoter state')
%grid on
    
%set(gca, 'Color', bkg_color)

promoter_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(promoter_fig,[figPath 'promoter_wt_trend.png'])
saveas(promoter_fig,[figPath 'promoter_wt_trend.pdf'])

initiation_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices)
    stairs(sweepResults(1).time_fine_cm/60,(sweepResults(1).promoter_state_array_cm(:,plot_indices(p))-1)*param_val_vec(7)*3/1e4...
            ,'Color',pm_color,'LineWidth',2)
end    
% plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
% scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 25])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('initiation rate (au/min)')
%grid on
    
%set(gca, 'Color', bkg_color)

initiation_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(initiation_fig,[figPath 'initiation_fig.png'])
saveas(initiation_fig,[figPath 'initiation_fig.pdf'])


ms2_fig = figure('Position',[100 100 512 256]);
hold on
for p = 1:length(plot_indices)
    plot(sweepInfo.time_axis_cm/60,sweepResults(1).ms2_traces_observed_cm(:,plot_indices(p))/1e4,'Color','k','LineWidth',3)
end    
% plot(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'Color','k','LineWidth',1)
% scatter(sweepInfo.time_axis_cm(1:3:end)/60,mean(sweepInfo.knirps_array_cm(1:3:end,:),2),'MarkerFaceColor',dark_green,'MarkerEdgeColor','k')
xlim([5 30])
ylim([0 45])
set(gca,'Fontsize',14)

xlabel('time (minutes)')
ylabel('promoter state')
%grid on
    
%set(gca, 'Color', bkg_color)

ms2_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% pbaspect([3 2 1])

saveas(ms2_fig,[figPath 'fluo_wt_trend.png'])
saveas(ms2_fig,[figPath 'fluo_wt_trend.pdf'])

%% Make a gallery of input-output plots
rng(567);
n_plot = 5;
n_options = size(sweepResults(1).knirps_traces_cm(:,plot_indices(p)),1);

close all

io_fig = figure;
tiledlayout(n_plot,n_plot);
plot_indices_2 = randsample(1:n_options,n_plot^2,false);

for i = 1:n_plot
    for  j = 1:n_plot
        ind = (i-1)*n_plot + j;
        nexttile;
        yyaxis left
        plot(sweepInfo.time_axis_cm/60,sweepResults(1).knirps_traces_cm(:,plot_indices_2(ind)),'Color',[dark_green 1],'LineWidth',2);
        set(gca,'ytick',[]);
        ylim([0 12])
        ax = gca;
        ax.YColor = 'k';

        yyaxis right
        noise_vec = normrnd(0,4,size(sweepInfo.time_axis_cm/60));
        plot(sweepInfo.time_axis_cm/60,sweepResults(1).ms2_traces_observed_cm(:,plot_indices_2(ind))/1e4+ noise_vec',...
                    '-k','LineWidth',1.5)
        set(gca,'xtick',[],'ytick',[]);
        xlim([5 30]);
        ylim([0 60])
        ax = gca;
        ax.YColor = 'k';

        box on
    end
end

saveas(io_fig,[figPath 'io_trace_panel.png'])
saveas(io_fig,[figPath 'io_trace_panel.pdf'])
