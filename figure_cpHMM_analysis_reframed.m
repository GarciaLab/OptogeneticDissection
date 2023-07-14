% script to generate figures summarizing cpHMM inference results
clear
close all

addpath(genpath('./lib'))

% Load data
%dateString = '03-Mar-2022 17-41-41';
%dateString = '02-Mar-2022 20-19-49';
%dateString = '12-Jun-2022 22-53-48'; % same number
%dateString = '12-Jun-2022 22-52-54'; % same number
%dateString = '12-Jun-2022 23-21-01'; % 4000
%dateString = '18-Jun-2022 12-32-48'; % 3500
%dateString = '18-Jun-2022 12-30-16'; % 3500
%dateString = '18-Jun-2022 12-28-04'; % 3000
dateString = '18-Jun-2022 12-18-19'; % 3000

%FolderString = 'HMM_result';
% dateString = '';
dataRoot = ['data' filesep 'burst_analysis_data' filesep];
%dataPath = [dataRoot FolderString filesep];
dataPath = [dataRoot dateString filesep];
%figPath = ['.' filesep 'fig' filesep 'burst_analysis_results_reframed' filesep FolderString filesep];
figPath = ['.' filesep 'fig' filesep 'burst_analysis_results_reframed' filesep dateString filesep];
mkdir(figPath)

load([dataPath filesep 'inference_summary.mat'],'inference_summary')

%% Generate basic figures summarizing parameter trends 
close all
knirps_offset = 0;
errorAlpha = 0.3; 

% set basic color parameters
cmap = brewermap(8,'Set2');
cmap_bu = brewermap(8,'Blues');
cmap_gr = brewermap(8,'Greens');
cmap_rd = brewermap(8,'Reds');
shape_cell = {'d','o','s'};
name_cell = {'low intensity', 'wildtype', 'high intensity'};
bkg_color = [228,221,209]/255;


% kon
kon_fig = figure;
hold on
s = [];
kon_vec_long = [];
kon_ste_vec_long = [];
knirps_vec_long = [];
knirps_ste_vec_long = [];
% kon_weight_vec_long = [];
  
for i = 1:length(inference_summary)
    % extract data
    knirps_vec = inference_summary(i).knirps_mean-knirps_offset;   
    kon_vec = inference_summary(i).kon_mean;
    nan_ft = ~isnan(kon_vec);
    kon_vec = kon_vec(nan_ft);
    knirps_vec = knirps_vec(nan_ft);
    knirps_vec_long = [knirps_vec_long knirps_vec];
    knirps_ste_vec_long = [knirps_ste_vec_long inference_summary(i).knirps_ste(nan_ft)]; 
    kon_vec_long = [kon_vec_long kon_vec];       
    kon_ste = inference_summary(i).kon_ste(nan_ft);
    kon_ste_vec_long = [kon_ste_vec_long kon_ste];
%     kon_weight_vec_long = [kon_weight_vec_long 1./kon_ste];
    
    % plot errorbars
    errorbar(knirps_vec,kon_vec,kon_ste,'.','Color',[0 0 0 .2],'Capsize',0)
    % scatter means
    s(end+1) = scatter(knirps_vec,kon_vec,50,shape_cell{i},'MarkerFaceColor',cmap_bu(i*2+1,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceAlpha',1);
end
% kon_weight_vec_long(isinf(kon_weight_vec_long)) = mean(kon_weight_vec_long(~isinf(kon_weight_vec_long)));
kon_weight_vec_long = 1./(kon_ste_vec_long + 0.01*nanmean(kon_vec_long));

set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
ylabel('burst frequency (events per minute)')
%grid on
    
legend(s,name_cell{:})
xlim([2 8])
ylim([0 3.8])
%set(gca, 'Color', bkg_color)

kon_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(kon_fig,[figPath 'kon_scatter.png'])
saveas(kon_fig,[figPath 'kon_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use bootstraping to conduct linear fit that weights each condition equally
% define weighted linear function to use for fitting
nBoots = 100;
knirps_axis = linspace(0,1.1*max(knirps_vec_long));
[kon_trend_mean, kon_trend_95, kon_trend_05, kon_trend_array] = fit_lin_trend(...
                                              knirps_axis,knirps_vec_long,knirps_ste_vec_long,kon_vec_long,...
                                                kon_ste_vec_long,kon_weight_vec_long.^2,nBoots);                                               

%{
%Jake: try binning
% plot the binned result
binNum_comp = 11;
binMin_comp = min(knirps_vec_long);
binMax_comp = max(knirps_vec_long);
edges_comp = linspace(binMin_comp,binMax_comp,binNum_comp);

[~,~,loc]=histcounts(knirps_vec_long,edges_comp);
meany = accumarray(loc(:),kon_vec_long(:))./accumarray(loc(:),1);
stdy = accumarray(loc(:),kon_vec_long(:),[],@std)./sqrt(accumarray(loc(:),1));
xmid = 0.5*(edges_comp(1:end-1)+edges_comp(2:end));

plot(xmid,meany,'-o')
%}


% add to plot
fill([knirps_axis fliplr(knirps_axis)],[kon_trend_95 fliplr(kon_trend_05)],cmap_bu(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p = plot(knirps_axis,kon_trend_mean,'Color',cmap_bu(7,:),'LineWidth',2);

legend([s p],[name_cell {'linear trend'}])

pbaspect([3 2 1])

saveas(kon_fig,[figPath 'kon_trend_scatter.png'])
saveas(kon_fig,[figPath 'kon_trend_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff

koff_fig = figure;
hold on
s = [];
koff_vec_long = [];
koff_ste_vec_long = [];
koff_weight_vec_long = [];
for i = 1:length(inference_summary)
    % extract data
    knirps_vec = inference_summary(i).knirps_mean-knirps_offset;      
    koff_vec = inference_summary(i).koff_mean;
    nan_ft = ~isnan(koff_vec);
    koff_vec = koff_vec(nan_ft);
    knirps_vec = knirps_vec(nan_ft);
    koff_vec_long = [koff_vec_long koff_vec];
    koff_ste = inference_summary(i).koff_ste(nan_ft);
    koff_ste_vec_long = [koff_ste_vec_long koff_ste];
%     koff_weight_vec_long = [koff_weight_vec_long 1./koff_ste];
    % plot errorbars
    errorbar(knirps_vec,koff_vec,koff_ste,'.','Color',[0 0 0 .2],'Capsize',0)
    % scatter means
    s(end+1) = scatter(knirps_vec,koff_vec,50,shape_cell{i},'MarkerFaceColor',cmap_gr(i*2+1,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceAlpha',1);
end
% koff_weight_vec_long(isinf(koff_weight_vec_long)) = mean(koff_weight_vec_long(~isinf(koff_weight_vec_long)));
koff_weight_vec_long = 1./(koff_ste_vec_long + 0.01*nanmean(koff_vec_long));

set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
ylabel('k_{off} (events per minute)')
%grid on
xlim([2 8])
legend(s,name_cell{:})

%set(gca, 'Color', bkg_color)

koff_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(koff_fig,[figPath 'koff_scatter.png'])
saveas(koff_fig,[figPath 'koff_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% use bootstraping to conduct linear fit that weights each condition equally
[koff_trend_mean, koff_trend_95, koff_trend_05, koff_trend_array] = fit_lin_trend(knirps_axis,knirps_vec_long,...
                                                    knirps_ste_vec_long,koff_vec_long,koff_ste_vec_long,koff_weight_vec_long.^2,nBoots);
% add to plot
fill([knirps_axis fliplr(knirps_axis)],[koff_trend_95 fliplr(koff_trend_05)],cmap_gr(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p = plot(knirps_axis,koff_trend_mean,'Color',cmap_gr(7,:),'LineWidth',2);

legend([s p],[name_cell {'linear trend'}])

pbaspect([3 2 1])

saveas(koff_fig,[figPath 'koff_trend_scatter.png'])
saveas(koff_fig,[figPath 'koff_trend_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%
% dur (1/koff)

dur_fig = figure;
hold on
s = [];
dur_vec_long = [];
dur_ste_vec_long = [];
% dur_weight_vec_long = [];
for i = 1:length(inference_summary)
    % extract data
    knirps_vec = inference_summary(i).knirps_mean-knirps_offset;    
    dur_vec = inference_summary(i).dur_mean/3600;
    nan_ft = ~isnan(dur_vec);
    dur_vec = dur_vec(nan_ft);
    knirps_vec = knirps_vec(nan_ft);
    
    dur_vec_long = [dur_vec_long dur_vec];          
    dur_ste = inference_summary(i).dur_ste(nan_ft)/3600;
    dur_ste_vec_long = [dur_ste_vec_long dur_ste];  
%     dur_weight_vec_long = [dur_weight_vec_long 1./dur_ste];
    
    % plot errorbars
    errorbar(knirps_vec,dur_vec,dur_ste,'.','Color',[0 0 0 .2],'Capsize',0)
    % scatter means
    s(end+1) = scatter(knirps_vec,dur_vec,50,shape_cell{i},'MarkerFaceColor',cmap_gr(i*2+1,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceAlpha',1);
end

% dur_weight_vec_long(isinf(dur_weight_vec_long)) = mean(dur_weight_vec_long(~isinf(dur_weight_vec_long)));
dur_weight_vec_long = 1./(dur_ste_vec_long + 0.01*nanmean(dur_vec_long));

set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
ylabel('burst duration (minutes)')
%grid on
xlim([2 8])
legend(s,name_cell{:},'Location','northwest')
ylim([0 0.9])
%set(gca, 'Color', bkg_color)

dur_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(dur_fig,[figPath 'dur_scatter.png'])
saveas(dur_fig,[figPath 'dur_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% use bootstraping to conduct linear fit that weights each condition equally
[dur_trend_mean, dur_trend_95, dur_trend_05, dur_trend_array] = fit_lin_trend(knirps_axis,knirps_vec_long,knirps_ste_vec_long,...
                                                              dur_vec_long,dur_ste_vec_long,dur_weight_vec_long.^2,nBoots);

% add to plot
fill([knirps_axis fliplr(knirps_axis)],[dur_trend_95 fliplr(dur_trend_05)],cmap_gr(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p = plot(knirps_axis,dur_trend_mean,'Color',cmap_gr(7,:),'LineWidth',2);

legend([s p],[name_cell {'linear trend'}])

pbaspect([3 2 1])

saveas(dur_fig,[figPath 'dur_trend_scatter.png'])
saveas(dur_fig,[figPath 'dur_trend_scatter.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%
% r
r_fig = figure;
hold on
s = [];
r_vec_long = [];
r_ste_vec_long = [];
r_weight_vec_long = [];

for i = 1:length(inference_summary)
    % extract data
    knirps_vec = inference_summary(i).knirps_mean-knirps_offset;
    r_vec = inference_summary(i).r2_mean*60;
    nan_ft = ~isnan(r_vec);
    r_vec = r_vec(nan_ft);
    knirps_vec = knirps_vec(nan_ft);
    r_vec_long = [r_vec_long r_vec];
    r_ste = inference_summary(i).r2_ste(nan_ft)*60;
    r_ste_vec_long = [r_ste_vec_long r_ste];
%     r_weight_vec_long = [r_weight_vec_long 1./r_ste];
    
    % plot errorbars
    errorbar(knirps_vec,r_vec,r_ste,'.','Color',[0 0 0 .2],'Capsize',0)
    % scatter means
    s(end+1) = scatter(knirps_vec,r_vec,50,shape_cell{i},'MarkerFaceColor',cmap_rd(i*2+1,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',1,'MarkerFaceAlpha',1);
end
% problem_flags = isinf(r_weight_vec_long);% | r_weight_vec_long <= prctile(r_weight_vec_long,5);
% r_weight_vec_long(problem_flags) = mean(~problem_flags);
r_weight_vec_long = 1./(r_ste_vec_long + 0.01*nanmean(r_vec_long));

set(gca,'Fontsize',14)

xlabel('[Knirps] (au)')
ylabel('burst amplitude (au per minute)')
%grid on
    
legend(s,name_cell{:},'Location','southwest')
xlim([2 8])
%set(gca, 'Color', bkg_color)
ylim([0 6e5])

r_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(r_fig,[figPath 'r_scatter.png'])
saveas(r_fig,[figPath 'r_scatter.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% use bootstraping to conduct linear fit that weights each condition equally
[r_trend_mean, r_trend_95, r_trend_05, r_trend_array] = fit_lin_trend(knirps_axis,knirps_vec_long,knirps_ste_vec_long,...
                                                              r_vec_long,r_ste_vec_long,r_weight_vec_long.^2,nBoots);

% add to plot
fill([knirps_axis fliplr(knirps_axis)],[r_trend_95 fliplr(r_trend_05)],cmap_rd(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p = plot(knirps_axis,r_trend_mean,'Color',cmap_rd(7,:),'LineWidth',2);

legend([s p],[name_cell {'linear trend'}])

pbaspect([3 2 1])

saveas(r_fig,[figPath 'r_trend_scatter.png'])
saveas(r_fig,[figPath 'r_trend_scatter.pdf'])

%% Construct a bar graph comparison

num_avg = 8;

[knirps_vec_long_sorted knirps_order] = sort(knirps_vec_long);

% kon
kon_avg_low = mean(kon_vec_long(knirps_order(1:num_avg)));
kon_std_low = std(kon_vec_long(knirps_order(1:num_avg)));
kon_avg_high = mean(kon_vec_long(knirps_order(end-num_avg+1:end)));
kon_std_high = std(kon_vec_long(knirps_order(end-num_avg+1:end)));

% duration
dur_avg_low = mean(dur_vec_long(knirps_order(1:num_avg)));
dur_std_low = std(dur_vec_long(knirps_order(1:num_avg)));
dur_avg_high = mean(dur_vec_long(knirps_order(end-num_avg+1:end)));
dur_std_high = std(dur_vec_long(knirps_order(end-num_avg+1:end)));

% koff
koff_avg_low = mean(koff_vec_long(knirps_order(1:num_avg)));
koff_std_low = std(koff_vec_long(knirps_order(1:num_avg)));
koff_avg_high = mean(koff_vec_long(knirps_order(end-num_avg+1:end)));
koff_std_high = std(koff_vec_long(knirps_order(end-num_avg+1:end)));

% r
r_avg_low = mean(r_vec_long(knirps_order(1:num_avg)));
r_std_low = std(r_vec_long(knirps_order(1:num_avg)));
r_avg_high = mean(r_vec_long(knirps_order(end-num_avg+1:end)));
r_std_high = std(r_vec_long(knirps_order(end-num_avg+1:end)));

% statistical test
%kon
x = kon_vec_long(knirps_order(1:num_avg));
y = kon_vec_long(knirps_order(end-num_avg+1:end));
[h p1] = ttest(x,y);

%dur
x = dur_vec_long(knirps_order(1:num_avg));
y = dur_vec_long(knirps_order(end-num_avg+1:end));
[h p2] = ttest(x,y);

fig = figure;
%tiledlayout(1,2)
%nexttile
bar([kon_avg_low/kon_avg_high 1;dur_avg_low/dur_avg_high 1; r_avg_high/r_avg_low 1])
ylim([0 2.75])
%bar([1 kon_avg_high/kon_avg_low;dur_avg_low/dur_avg_high 1; 1 r_avg_high/r_avg_low])
%bar([kon_avg_low kon_avg_high;dur_avg_low dur_avg_high])
%bar([kon_avg_low kon_avg_high])
%bar([kon_avg_low kon_avg_high;koff_avg_low koff_avg_high])
%bar([1 kon_avg_high/kon_avg_low;1 koff_avg_high/koff_avg_low; 1 r_avg_high/r_avg_low])
%nexttile
%bar([1 kon_avg_high/kon_avg_low;1 dur_avg_high/dur_avg_low; 1 r_avg_high/r_avg_low])

fig = figure;
tiledlayout(1,3)
%daspect([3 1 1]);
%pbaspect([1 2 1]);
nexttile
hold on
bar([kon_avg_low kon_avg_high])
errorbar([kon_avg_low kon_avg_high],[kon_std_low kon_std_high ])
ylim([0 2.5*kon_avg_high])
box on


nexttile
hold on
bar([dur_avg_low dur_avg_high])
errorbar([dur_avg_low dur_avg_high],[dur_std_low dur_std_high ])
%ylim([0 2.6*dur_avg_low])
ylim([0 2*dur_avg_high])
box on
%set(gca, 'YScale', 'log')


nexttile
hold on
bar([r_avg_low r_avg_high])
errorbar([r_avg_low r_avg_high],[r_std_low r_std_high], '. ')
ylim([0 2*r_avg_high])
box on
%xlim([0 5])
% try log plot
%set(gca, 'YScale', 'log')

%pbaspect([3 2 1])


%% Attempt to generate summary figures
close all
% Full figure with everything

% normalize predictions by values at Knirps = 3AU
[~,norm_ind] = min(abs(knirps_axis-3));
kon_norm = kon_trend_mean(norm_ind);
koff_norm = koff_trend_mean(norm_ind);
dur_norm = dur_trend_mean(norm_ind);
r_norm = r_trend_mean(norm_ind);

markerAlpha = 0.25;
errorAlpha = 0.5;
markerSize = 40;
full_figure = figure;
hold on

for i = 1:length(inference_summary)
    % extract data
    knirps_vec = inference_summary(i).knirps_mean-knirps_offset;    
    dur_vec = inference_summary(i).dur_mean / dur_norm / 3600;
    kon_vec = inference_summary(i).kon_mean / kon_norm;
    koff_vec = inference_summary(i).koff_mean / koff_norm;
    r_vec = inference_summary(i).r2_mean * 60 / r_norm;    
        
    % get error
    r_ste = inference_summary(i).r2_ste*60 / r_norm;
    kon_ste = inference_summary(i).kon_ste / kon_norm;
    koff_ste = inference_summary(i).koff_ste / koff_norm;
    dur_ste = inference_summary(i).dur_ste / dur_norm / 3600;
    
    % plot errorbars
%     e1 = errorbar(knirps_vec,r_vec,r_ste,'.','MarkerSize',0.1,'Color','k','Capsize',0);
%     set([e1.Bar, e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*0.1])
%     e1 = errorbar(knirps_vec,dur_vec,dur_ste,'.','MarkerSize',0.1,'Color','k','Capsize',0);
%     set([e1.Bar, e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*0.1])
%     e1 = errorbar(knirps_vec,kon_vec,kon_ste,'.','MarkerSize',0.1,'Color','k','Capsize',0);
%     set([e1.Bar, e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*0.1])
    
    % scatter means
    scatter(knirps_vec,kon_vec,markerSize,shape_cell{1},'MarkerFaceColor',cmap_bu(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',markerAlpha,'MarkerFaceAlpha',markerAlpha);
    %scatter(knirps_vec,koff_vec,markerSize,shape_cell{1},'MarkerFaceColor',cmap_gr(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',markerAlpha,'MarkerFaceAlpha',markerAlpha);
    scatter(knirps_vec,dur_vec,markerSize,shape_cell{1},'MarkerFaceColor',cmap_gr(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',markerAlpha,'MarkerFaceAlpha',markerAlpha);
    scatter(knirps_vec,r_vec,markerSize,shape_cell{1},'MarkerFaceColor',cmap_rd(5,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',markerAlpha,'MarkerFaceAlpha',markerAlpha);
end

% plot linear trends 
fill([knirps_axis fliplr(knirps_axis)],[kon_trend_95 fliplr(kon_trend_05)]/kon_norm,cmap_bu(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p1 = plot(knirps_axis,kon_trend_mean/kon_norm,'Color',cmap_bu(7,:),'LineWidth',2);

%fill([knirps_axis fliplr(knirps_axis)],[koff_trend_95 fliplr(koff_trend_05)]/koff_norm,cmap_gr(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
%p2 = plot(knirps_axis,_trend_mean/koff_norm,'Color',cmap_gr(7,:),'LineWidth',2);

fill([knirps_axis fliplr(knirps_axis)],[dur_trend_95 fliplr(dur_trend_05)]/dur_norm,cmap_gr(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p2 = plot(knirps_axis,dur_trend_mean/dur_norm,'Color',cmap_gr(7,:),'LineWidth',2);

fill([knirps_axis fliplr(knirps_axis)],[r_trend_95 fliplr(r_trend_05)]/r_norm,cmap_rd(2,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p3 = plot(knirps_axis,r_trend_mean/r_norm,'Color',cmap_rd(7,:),'LineWidth',2);

set(gca,'Fontsize',14)

% try log plot
set(gca, 'YScale', 'log')
%ylim([0 4])

xlabel('[Knirps] (au)')
ylabel('normalized burst parameter trends')
%grid on

xlim([2 9])
legend([p1 p2 p3],'burst frequency (k_{on})','burst duration (1/k_{off})','burst amplitude (r)','Location','northwest')
ylim([0.1 10])
%set(gca, 'Color', bkg_color)

full_figure.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(full_figure,[figPath 'all_trend_scatter.png'])
saveas(full_figure,[figPath 'all_trend_scatter.pdf'])

%% Plot relative contributions to repression
errorAlpha = 0.3;
cmap_gra = brewermap(8,'Set2');

% calculate trend for raw data
data_trend = r_vec_long.*kon_vec_long./(kon_vec_long + koff_vec_long);

% full trend
predicted_rate_array = r_trend_array .* kon_trend_array ./ (kon_trend_array + koff_trend_array);
predicted_rate_95 = prctile(predicted_rate_array,95,2)';
predicted_rate_05 = prctile(predicted_rate_array,5,2)';

predicted_rate_trend = r_trend_mean.*kon_trend_mean./(kon_trend_mean + koff_trend_mean);

% no kon modulation
predicted_rate_array_no_kon = r_trend_array .* kon_trend_array(norm_ind,:) ./ (kon_trend_array(norm_ind,:) + koff_trend_array);
predicted_rate_no_kon_95 = prctile(predicted_rate_array_no_kon,95,2)';
predicted_rate_no_kon_05 = prctile(predicted_rate_array_no_kon,5,2)';

rate_prediction_no_kon = r_trend_mean.*kon_trend_mean(norm_ind)./(kon_trend_mean(norm_ind) + koff_trend_mean);

% only kon modulation
predicted_rate_array_kon_only = r_trend_array(norm_ind,:) .* kon_trend_array ./ (kon_trend_array + koff_trend_array(norm_ind,:));
predicted_rate_kon_only_95 = prctile(predicted_rate_array_kon_only,95,2)';
predicted_rate_kon_only_05 = prctile(predicted_rate_array_kon_only,5,2)';

rate_prediction_kon_only = r_trend_mean(norm_ind).*kon_trend_mean./(kon_trend_mean + koff_trend_mean(norm_ind));


repression_figure = figure;
hold on

plot(knirps_axis,repelem(predicted_rate_trend(norm_ind),length(knirps_axis)),':','Color','k','LineWidth',2)

% fill([knirps_axis fliplr(knirps_axis)],[predicted_rate_no_kon_95 fliplr(predicted_rate_no_kon_05)],cmap_gra(3,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p2 = plot(knirps_axis,rate_prediction_no_kon,'--','Color',cmap_gra(2,:),'LineWidth',2);

% fill([knirps_axis fliplr(knirps_axis)],[predicted_rate_kon_only_95 fliplr(predicted_rate_kon_only_05)],cmap_gra(8,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p3 = plot(knirps_axis,rate_prediction_kon_only,'-.','Color',cmap_gra(8,:),'LineWidth',2);

fill([knirps_axis fliplr(knirps_axis)],[predicted_rate_95 fliplr(predicted_rate_05)],cmap_gra(3,:),'FaceAlpha',errorAlpha,'EdgeAlpha',errorAlpha)
p1 = plot(knirps_axis,predicted_rate_trend,'Color',cmap_gra(3,:),'LineWidth',2);

scatter(knirps_vec_long,data_trend,50,'d','MarkerFaceColor',cmap_gra(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5)

set(gca,'Fontsize',14)

% try log plot
%set(gca, 'YScale', 'log')

xlabel('[Knirps] (au)')
ylabel({'predicted transcription rate';'(au per minute)'})
%grid on

xlim([2 8])
legend([p1 p2 p3],'full model','no frequency modulation','only frequency modulation','Location','southwest')
% ylim([0 2])
%set(gca, 'Color', bkg_color)

repression_figure.InvertHardcopy = 'off';
set(gcf,'color','w');

pbaspect([3 2 1])

saveas(repression_figure,[figPath 'rate_trend_scatter.png'])
saveas(repression_figure,[figPath 'rate_trend_scatter.pdf'])