% script to generate figures summarizing cpHMM inference results
clear
close all

addpath(genpath('./lib'))

% Load data
dateString = '28-Jun-2022 15-05-19'; 
dataRoot = ['data' filesep 'burst_analysis_data' filesep];
dataPath = [dataRoot dateString filesep];

figPath = ['.' filesep 'fig' filesep 'burst_analysis_results_binary' filesep dateString filesep];
mkdir(figPath)

load([dataPath filesep 'inference_summary.mat'],'inference_summary')

%% use bootstrap sampling to construct delta distributions 


kni_group_vec = inference_summary(1).knirps_group;
kon_vec_long = inference_summary(1).kon_vec;
dur_vec_long = inference_summary(1).dur_vec;
r2_vec_long = inference_summary(1).r2_vec;
[h_kon,p_kon] = kstest2(kon_vec_long(kni_group_vec==1),kon_vec_long(kni_group_vec==3))
[h_dur,p_dur] = kstest2(dur_vec_long(kni_group_vec==1),dur_vec_long(kni_group_vec==3))
[h_r,p_r] = ttest2(r2_vec_long(kni_group_vec==1),r2_vec_long(kni_group_vec==3),'Vartype','unequal')
%%
rng(500)
nBoots = 1e4;
kon_delta_vec = NaN(1,nBoots);
dur_delta_vec = NaN(1,nBoots);
r_delta_vec = NaN(1,nBoots);
k1_options = find(kni_group_vec==1);
k3_options = find(kni_group_vec==3);
for n = 1:nBoots
    k1_samps = randsample(k1_options,1);
    k3_samps = randsample(k3_options,1);
    kon_delta_vec(n) = kon_vec_long(k1_samps) - kon_vec_long(k3_samps);
    dur_delta_vec(n) = mean(dur_vec_long(k3_samps)) - dur_vec_long(k1_samps);
    r_delta_vec(n) = r2_vec_long(k1_samps) - r2_vec_long(k3_samps);
end  

kon_delta_mean = mean(kon_delta_vec);
kon_delta_std = std(kon_delta_vec);
kon_p_boot = mean(kon_delta_vec<=0)
kon_p_alt = normcdf(0,kon_delta_mean,kon_delta_std)

dur_delta_mean = mean(dur_delta_vec);
dur_delta_std = std(dur_delta_vec);
dur_p_boot = mean(dur_delta_vec<=0) % note that duration p hovers just beyond 10% significance cut-off
dur_p_alt = normcdf(0,dur_delta_mean,dur_delta_std)

r_delta_mean = mean(r_delta_vec);
r_delta_std = std(r_delta_vec);
r_p_boot = mean(r_delta_vec<=0)
r_p_alt = normcdf(0,r_delta_mean,r_delta_std)
%% Make high-vs.-low bar plots for burst parameters

cmap = brewermap([],'Set2');
close all

kon_vec_mean = inference_summary(1).kon_mean;
kon_vec_ste = inference_summary(1).kon_ste;

dur_vec_mean = inference_summary(1).dur_mean;
dur_vec_ste = inference_summary(1).dur_ste;

r_vec_mean = inference_summary(1).r2_mean;
r_vec_ste = inference_summary(1).r2_ste;

burst_fig = figure('Position',[100 100 512 256]);
tiledlayout(1,3)
nexttile
hold on
bar(1,kon_vec_mean(1),'FaceColor',cmap(7,:))
bar(2,kon_vec_mean(3),'FaceColor',cmap(8,:))
errorbar(1:2,kon_vec_mean([1 3]),kon_vec_ste([1 3]),'.','Color','k','LineWidth',1.25)
ylim([0 1.5*max(kon_vec_mean)])
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on

nexttile
hold on
bar(1,dur_vec_mean(1),'FaceColor',cmap(7,:))
bar(2,dur_vec_mean(3),'FaceColor',cmap(8,:))
errorbar(1:2,dur_vec_mean([1 3]),dur_vec_ste([1 3]),'.','Color','k','LineWidth',1.25)
%ylim([0 2.6*dur_avg_low])
ylim([0 1.5*max(dur_vec_mean)])
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on

nexttile
hold on
bar(1,r_vec_mean(1),'FaceColor',cmap(7,:))
bar(2,r_vec_mean(3),'FaceColor',cmap(8,:))
errorbar(1:2,r_vec_mean([1 3]),r_vec_ste([1 3]),'.','Color','k','LineWidth',1.25)
ylim([0 1.5*max(r_vec_mean)])
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on
saveas(burst_fig,[figPath 'burst_bar_binary.pdf'])