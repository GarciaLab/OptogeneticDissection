clear all
close all
clc

%% create color map

%Magenta color map
cmap_magenta = [217 82 147];
cmap_white = [255 255 255];

new_stepNum = 256;

cmap_magenta_1 = interp1([0 1],[double(cmap_white(1)) double(cmap_magenta(1))],linspace(0,1,new_stepNum));
cmap_magenta_2 = interp1([0 1],[double(cmap_white(2)) double(cmap_magenta(2))],linspace(0,1,new_stepNum));
cmap_magenta_3 = interp1([0 1],[double(cmap_white(3)) double(cmap_magenta(3))],linspace(0,1,new_stepNum));

cmap_magenta_new = [cmap_magenta_1' cmap_magenta_2' cmap_magenta_3']/256;

beta = -0.5;
cmap_magenta_updated = brighten(cmap_magenta_new,beta);

%%

wildtype = load('data/eve_rescue/2020-01-16-eveMS2_control_embryo1/CompiledParticles.mat');
optokni = load('data/eve_rescue/2019-12-12-eYFP_homo_embryo3/CompiledParticles.mat');

wildtype_nc14 = 20;
optokni_nc14 = 1;

wildtype.ElapsedTime = wildtype.ElapsedTime-wildtype.ElapsedTime(wildtype_nc14);
optokni.ElapsedTime = optokni.ElapsedTime-optokni.ElapsedTime(optokni_nc14);

%% Interpolate data for comparison

timeRange = linspace(0,40,81);

wildtype.MeanRate = wildtype.MeanVectorAP{1,1}.*wildtype.OnRatioAP{1,1};
optokni.MeanRate = optokni.MeanVectorAP{1,1}.*optokni.OnRatioAP{1,1};

wildtype.OnRatio = wildtype.OnRatioAP{1,1};
optokni.OnRatio = optokni.OnRatioAP{1,1};

wildtype.MeanRateInterp = interp1(wildtype.ElapsedTime,wildtype.MeanRate,timeRange);
optokni.MeanRateInterp = interp1(optokni.ElapsedTime,optokni.MeanRate,timeRange);

wildtype.OnRatioInterp = interp1(wildtype.ElapsedTime,wildtype.OnRatio,timeRange);
optokni.OnRatioInterp = interp1(optokni.ElapsedTime,optokni.OnRatio,timeRange);

%% Plot interpolated data

wildtype.APbinID = wildtype.APbinID - 0.050;


fig = figure;
imagesc(wildtype.APbinID, timeRange, wildtype.OnRatioInterp)
%colormap(cmap_magenta_updated)
colorbar
xlim([0.47 0.72])
ylim([0 40])
clim([0 1])
pbaspect([3 2 1])
xlabel('anterior-posterior position (% embryo length)')
ylabel('time (min)')

fig = figure;
imagesc(optokni.APbinID, timeRange, optokni.OnRatioInterp)
%colormap(cmap_magenta_updated)
colorbar
xlim([0.47 0.72])
ylim([0 40])
clim([0 1])
pbaspect([3 2 1])
xlabel('anterior-posterior position (% embryo length)')
ylabel('time (min)')


%% plot transcription rate

time_plot = 75;

fig = figure;
hold on
plot(wildtype.APbinID,wildtype.OnRatioInterp(time_plot,:))
plot(optokni.APbinID,optokni.OnRatioInterp(time_plot,:))

xlim([0.47 0.72])

fig = figure;
hold on
plot(wildtype.APbinID,wildtype.MeanRateInterp(time_plot,:))
plot(optokni.APbinID,optokni.MeanRateInterp(time_plot,:))

xlim([0.47 0.72])

%% calculate integrated data

k = 7; % mRNA half life
dt = timeRange(2)-timeRange(1);

wildtype.mRNA_int = zeros(length(timeRange),length(wildtype.APbinID));
optokni.mRNA_int = zeros(length(timeRange),length(optokni.APbinID));

for i = 2:length(timeRange)

    MeanRate_wt = wildtype.MeanRateInterp(i,:);
    MeanRate_wt(isnan(MeanRate_wt)) = 0;

    MeanRate_opto = optokni.MeanRateInterp(i,:);
    MeanRate_opto(isnan(MeanRate_opto)) = 0;
    
    % with degradation
    wildtype.mRNA_int(i,:) = wildtype.mRNA_int(i-1,:)*2^(-dt/k) + MeanRate_wt*dt;
    optokni.mRNA_int(i,:) = optokni.mRNA_int(i-1,:)*2^(-dt/k) + MeanRate_opto*dt;
    
end

%% Plot integrated mRNA

fig = figure;
imagesc(wildtype.APbinID,timeRange,wildtype.mRNA_int)

fig = figure;
imagesc(optokni.APbinID,timeRange,optokni.mRNA_int)

norm_factor = 1.6659;
amp_max = 393392;

time_plot = [41 61 81];

for i = 1:3
    fig = figure;
    hold on
    plot(wildtype.APbinID,movmean(wildtype.mRNA_int(time_plot(i),:),3)/amp_max,'- .','LineWidth',2,'MarkerSize',15)
    plot(optokni.APbinID,movmean(optokni.mRNA_int(time_plot(i),:)*norm_factor/amp_max,3),'- .','LineWidth',2,'MarkerSize',15)
    xlabel('anterior-posterior position (% embryo length)')
    ylabel('integrated mRNA (normalized)')
    %xlim([0.47 0.72])
    ylim([0 1.2])
    pbaspect([3 2 1])
end