clear
close all
clc

addpath(genpath('./lib'))

%% initialization

projectName = 'optokni_eve4+6_ON'; 

% load data
load(['data' filesep 'main_analysis' filesep projectName filesep 'spot_struct.mat'])

% specify the parameters for each embryo
embryo(1).expID = 1;
embryo(1).frame_on = 36;

embryo(2).expID = 2;
embryo(2).frame_on = 43;

embryo(3).expID = 3;
embryo(3).frame_on = 45;

embryo(4).expID = 4;
embryo(4).frame_on = 43;

% color to be used
%color_green = brighten([38 142 75]/256,.4);
color_green = [122 168 116]/255;
mRNA_red = brighten([212 100 39]/256,.2);
magenta = '#D95292';

knirps_offset = 375698.13;

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

time_threshold = 2; %min

% histogram parameters
binNum = 15;
binMax = 8;
edges = linspace(0,binMax,binNum);

% timerange to analyze for response time
analysis_range = 8;


%% Figure: plot mean fluorescence vs time

data_filter_full = [];
silence_time_full = [];
response_time_full = []; 

time_aligned_full_long = [];
fluo_full_long = [];
knirps_full_long = [];
frame_full_long = [];
first_on_all_time_long = [];
last_on_all_time_long = [];

temp_traj_fig  = figure('Position',[0 0 500 500]);

nuclei_count = 0;

for i = 1:length(embryo)

    expID = embryo(i).expID;
    frame_on = embryo(i).frame_on;
    
    nBoots = 100;

    ever_on_vec = [];
    mean_ap = [];
    time_orig_long = [];
    fluo_orig_long = [];
    frame_orig_long = [];
    off_time_long = [];
    knirps_orig_long = [];
    last_on_long = [];
    first_on_long = [];
    first_on_all_long = [];
    last_on_all_long = [];

    for j = 1:length(spot_struct)

        if (spot_struct(j).TraceQCFlag == 1) && (spot_struct(j).setID == expID)
            % extract core vectors 

            % extract core vectors 
            fluo_vec_orig = spot_struct(j).fluo;
            time_vec_orig = spot_struct(j).time;
            frame_vec_orig = spot_struct(j).frames;
            knirps_vec_orig = spot_struct(j).rawNCProtein;
            ap_vec_orig = spot_struct(j).APPosNucleus;

            % calculate mean
            ever_on_orig = any(~isnan(fluo_vec_orig));
            mean_ap_orig = nanmean(ap_vec_orig);
            mean_knirps_orig = nanmean(knirps_vec_orig);

            if ever_on_orig
                last_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig <= frame_on),1,'last'));
                first_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig > frame_on),1));
                
            end


            if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)

               nuclei_count = nuclei_count + 1;
               
               time_orig_long = [time_orig_long time_vec_orig];
               frame_orig_long = [frame_orig_long frame_vec_orig];
               fluo_orig_long = [fluo_orig_long fluo_vec_orig];
               knirps_orig_long = [knirps_orig_long knirps_vec_orig];

               if isempty(first_on_frame)
                   first_on_all_long = [first_on_all_long NaN];
                   last_on_all_long = [last_on_all_long NaN];
               end

               if ~isempty(last_on_frame) && ~isempty(first_on_frame)
                   last_on_long = [last_on_long last_on_frame];
                   last_on_all_long = [last_on_all_long last_on_frame];
                   first_on_long = [first_on_long first_on_frame];
                   first_on_all_long = [first_on_all_long first_on_frame];
               end
            end

        end

    end
    
    
    % calculate mean knirps and fraction on (before/after perturbation)

    fluo_orig_long(isnan(fluo_orig_long)) = 0;

    frame_len = max(frame_orig_long);
    
    time_vec = zeros(1,frame_len);
    fluo_vec_mean = zeros(frame_len,1);
    fluo_vec_ste = zeros(frame_len,1);
    knirps_vec_mean = zeros(frame_len,1);
    knirps_vec_ste = zeros(frame_len,1);
    frac_on = zeros(frame_len,1);
    frac_on_ste = zeros(frame_len,1);

    fluo_orig_long_zero = fluo_orig_long;
    fluo_orig_long_zero(isnan(fluo_orig_long)) = 0;

    fluo_orig_long_binary = fluo_orig_long;
    fluo_orig_long_binary(isnan(fluo_orig_long)) = 0;
    fluo_orig_long_binary(fluo_orig_long>0) = 1;


    for j = 1:frame_len

        time_filter_long = frame_orig_long==j;
        time_vec(j) = mean(time_orig_long(time_filter_long))/60;

        boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long_zero(time_filter_long));
        fluo_vec_mean(j) = nanmean(boot_samples_fluo);
        fluo_vec_ste(j) = std(boot_samples_fluo,'omitnan');
 
        boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_orig_long(time_filter_long));
        knirps_vec_mean(j) = nanmean(boot_samples_knirps);
        knirps_vec_ste(j) = std(boot_samples_knirps,'omitnan');

        frac_on(j) = nanmean(fluo_orig_long_binary(time_filter_long));
        frac_on_ste(j) = std(fluo_orig_long_binary(time_filter_long));

    end

    time_vec_on = time_vec-time_vec(frame_on);

    %knirps_vec_mean(time_vec_on>=0) = knirps_vec_mean(time_vec_on>=0)/correction_factor;
    knirps_vec_mean(time_vec_on>=0) = calibrate_export_laser_custom(knirps_vec_mean(time_vec_on>=0));    
    knirps_vec_mean = knirps_vec_mean-knirps_offset;

    % record the result for this embryo
    data_filter = (time_vec(last_on_long) <= time_vec(frame_on)-time_threshold);
    data_filter_full = [data_filter_full double(data_filter)];
    
    silence_time_full = [silence_time_full time_vec(frame_on)-time_vec(last_on_long)];
    response_time_full = [response_time_full time_vec(first_on_long)-time_vec(frame_on)];

    for j = 1:length(first_on_all_long)
        if first_on_all_long(j) == 0
            first_on_all_time_long = [first_on_all_time_long NaN];
        else
            first_on_all_time_long = [first_on_all_time_long time_vec_on(first_on_all_long(j))];
        end
    end

    for j = 1:length(last_on_all_long)
        if last_on_all_long(j) == 0
            last_on_all_time_long = [last_on_all_time_long NaN];
        else
            last_on_all_time_long = [last_on_all_time_long time_vec_on(last_on_all_long(j))];
        end
    end
    
    % record results for combining the traces
    time_vec_aligned = time_orig_long/60 - time_vec(frame_on);
    time_aligned_full_long = [time_aligned_full_long time_vec_aligned];
    fluo_full_long = [fluo_full_long fluo_orig_long];
    knirps_full_long = [knirps_full_long knirps_orig_long];
    frame_full_long = [frame_full_long frame_orig_long];
    
    
    hold on
    time_interp = -10:0.1:10;
    frac_on_mean = movmean(frac_on,5);
    frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
    frac_on_interp = movmean(frac_on_interp,5);

    yyaxis left

    errorbar(time_vec_on,knirps_vec_mean,knirps_vec_ste,'Color','k','CapSize',0);
    plot(time_vec_on,knirps_vec_mean,'-k','LineWidth',1)
    switch i
        case 1
            scatter(time_vec_on,knirps_vec_mean,75,'MarkerFaceColor',color_green,'MarkerEdgeColor','k')
        case 2
            scatter(time_vec_on,knirps_vec_mean,90,'s','MarkerFaceColor',color_green,'MarkerEdgeColor','k')
        case 3
            scatter(time_vec_on,knirps_vec_mean,75,'d','MarkerFaceColor',color_green,'MarkerEdgeColor','k')
        case 4
            scatter(time_vec_on,knirps_vec_mean,75,'^','MarkerFaceColor',color_green,'MarkerEdgeColor','k')
    end

    ylim([3.75E5 9E5])
    ylabel(['Knirps concentration (AU)'])

    yyaxis right
    hold on
    plot(time_vec_on,frac_on,'.')
    plot(time_interp,frac_on_interp,'-','LineWidth',2,'Color',magenta);
    switch i
        case 1
            scatter(time_vec_on,frac_on,75,'MarkerFaceColor',magenta,'MarkerEdgeColor','k')
        case 2
            scatter(time_vec_on,frac_on,90,'s','MarkerFaceColor',magenta,'MarkerEdgeColor','k')
        case 3
            scatter(time_vec_on,frac_on,75,'d','MarkerFaceColor',magenta,'MarkerEdgeColor','k')
        case 4
            scatter(time_vec_on,frac_on,75,'^','MarkerFaceColor',magenta,'MarkerEdgeColor','k')
    end
    ylim([0 1])
    ylabel(['fraction of nuclei on'])

    ax = gca;
    ax.YAxis(1).Color = color_green;
    ax.YAxis(2).Color = magenta;

    xlim([-10 7])
    xlabel(['time relative to perturbation (min)'])

    pbaspect([3 2 1])

end

%% plot combined results

% calculate mean knirps and fraction on (before/after perturbation)

fluo_full_long(isnan(fluo_full_long)) = 0;

time_bin_full = linspace(-15,15,58);
%time_bin_full = linspace(-15,15,81);

time_vec_plot = (time_bin_full(1:end-1) + time_bin_full(2:end))/2;
time_groups_full = discretize(time_aligned_full_long,time_bin_full);

fluo_vec_full_mean = zeros(length(time_bin_full)-1,1);
fluo_vec_full_ste = zeros(length(time_bin_full)-1,1);

knirps_vec_full_mean = zeros(length(time_bin_full)-1,1);
knirps_vec_full_ste = zeros(length(time_bin_full)-1,1);

frac_on_full = zeros(length(time_bin_full)-1,1);
frac_on_full_ste = zeros(length(time_bin_full)-1,1);

fluo_full_long_zero = fluo_full_long;
fluo_full_long_zero(isnan(fluo_full_long)) = 0;

fluo_full_long_binary = fluo_full_long;
fluo_full_long_binary(isnan(fluo_full_long)) = 0;
fluo_full_long_binary(fluo_full_long>0) = 1;


for j = 1:length(time_bin_full)-1

    time_filter_long_full = time_groups_full==j;

    time_vec_full(j) = mean(time_aligned_full_long(time_filter_long_full));

    boot_samples_fluo_full = bootstrp(nBoots,@nanmean,fluo_full_long_zero(time_filter_long_full));
    fluo_vec_full_mean(j) = nanmean(boot_samples_fluo_full);
    fluo_vec_full_ste(j) = std(boot_samples_fluo_full,'omitnan');
 
    boot_samples_knirps_full = bootstrp(nBoots,@nanmean,knirps_full_long(time_filter_long_full));
    knirps_vec_full_mean(j) = nanmean(boot_samples_knirps_full);
    knirps_vec_full_ste(j) = std(boot_samples_knirps_full,'omitnan');

    boot_samples_frac_on_full = bootstrp(nBoots,@nanmean,fluo_full_long_binary(time_filter_long_full));
    frac_on_full(j) = nanmean(fluo_full_long_binary(time_filter_long_full));
    %frac_on_full_ste(j) = std(fluo_full_long_binary(time_filter_long_full));
    frac_on_full_ste(j) = std(boot_samples_frac_on_full,'omitnan');


end

fig = figure;

knirps_vec_full_mean(time_vec_plot>=0) = calibrate_export_laser_custom(knirps_vec_full_mean(time_vec_plot>=0));
knirps_vec_full_mean = knirps_vec_full_mean - knirps_offset;


hold on
time_full_interp = -10:0.1:10;
frac_on_full_mean = movmean(frac_on_full,3);
frac_on_full_interp = interp1(time_vec_full(~isnan(frac_on_full_mean)),frac_on_full_mean(~isnan(frac_on_full_mean)),time_full_interp,'spline');
frac_on_full_interp = movmean(frac_on_full_interp,1);

yyaxis left
errorbar(time_vec_plot,knirps_vec_full_mean,knirps_vec_full_ste,'Color','k','CapSize',0);
plot(time_vec_plot,knirps_vec_full_mean,'-k','LineWidth',1)
scatter(time_vec_plot,knirps_vec_full_mean,90,'MarkerFaceColor',color_green,'MarkerEdgeColor','k')
ylim([3.75E5 9E5])
ylabel(['Knirps concentration (AU)'])


yyaxis right
%plot(time_vec_plot,frac_on_full,'-','LineWidth',1,'Color','k')
errorbar(time_vec_plot,frac_on_full,frac_on_full_ste,'.','Color','k','CapSize',0)
plot(time_full_interp,frac_on_full_interp,'-','LineWidth',2,'Color',magenta);
scatter(time_vec_plot,frac_on_full,90,'MarkerFaceColor',magenta,'MarkerEdgeColor','k')
ylabel(['fraction of nuclei on'])
ylim([0 0.9])

ax = gca;
ax.YAxis(1).Color = color_green;
ax.YAxis(2).Color = magenta;


xlim([-10 7])
xlabel(['time relative to perturbation (min)'])

pbaspect([3 2 1])

%% plot single traces

sample_traces = [];
time_vec_interp = -10:0.1:7;

% extract single-trace data
for j = 1:length(spot_struct)
    
    embryo_num = find(spot_struct(j).setID == [embryo(:).expID]);
    
    if ~isnan(embryo_num)
        frame_on = embryo(embryo_num).frame_on;

        ap_vec = spot_struct(j).APPosNucleus;
        ap_pos = mean(ap_vec);

        if (ap_pos>=-0.02) && (ap_pos<=0.02)
            spot_fluo = spot_struct(j).fluo;
            knirps_fluo = spot_struct(j).rawNCProtein;
            spot_fluo(isnan(spot_fluo)) = 0;
            
            time_vec = spot_struct(j).time/60;
            frame_vec = spot_struct(j).frames;

            frame_start = frame_vec(1);
            frame_final = frame_vec(end);

            if (frame_on>=frame_start) && (frame_on<=frame_final) && ~isnan(frame_start) && ~isnan(frame_final)
                time_on = time_vec(frame_vec == frame_on);
                time_vec_new = time_vec-time_on;

                spot_fluo_interp = interp1(time_vec_new,spot_fluo,time_vec_interp,'linear');

                sample_traces = [sample_traces;spot_fluo_interp]; 

            end
        end
    end
    
end

last_on_time_long = zeros(size(sample_traces,1),1);
first_on_time_long = zeros(size(sample_traces,1),1);

for i = 1:size(sample_traces,1)
    spot_vec_temp = sample_traces(i,:);
    
    last_on_time = time_vec_interp(find((spot_vec_temp>0) & (time_vec_interp<0),1,'last'));
    first_on_time = time_vec_interp(find((spot_vec_temp>0) & (time_vec_interp>0),1));
    
    if ~isempty(last_on_time) && ~isempty(first_on_time) && (last_on_time<-2)
        last_on_time_long(i) = last_on_time;
        first_on_time_long(i) = first_on_time;
    else
        last_on_time_long(i) = NaN;
        first_on_time_long(i) = NaN;
    end
end

[B,I] = sort(first_on_time_long,'descend');

%I = I((~isnan(B)) & (B<=5));
I = I((~isnan(B)));

sample_traces_fig = figure;
imagesc('XData',time_vec_interp,'CData',sample_traces(I,:)*1e-5)
xlim([-10 7])
ylim([1 length(I)])
xlabel('time relative to perturbation (min)')
ylabel('trace num')
colormap(plasma)
caxis([0 4.5])
colorbar
pbaspect([3 2 1])

%% calculate silenced duration vs response time

analysis_range_sil_dur = 12;

response_time_final = response_time_full(data_filter_full==1);
silence_time_final = silence_time_full(data_filter_full==1);


x = silence_time_final;
y = response_time_final;
filter = (x<analysis_range_sil_dur) & (y<analysis_range);

xFit = x(filter);
yFit = y(filter);


% plot the binned result
binNum_comp = 9;
binMax_comp = 11.5;
edges_comp = linspace(time_threshold,binMax_comp,binNum_comp);

[~,~,loc]=histcounts(xFit,edges_comp);
meany = accumarray(loc(:),yFit(:))./accumarray(loc(:),1);
stdy = accumarray(loc(:),yFit(:),[],@std)./sqrt(accumarray(loc(:),1));
xmid = 0.5*(edges_comp(1:end-1)+edges_comp(2:end));

test_fig = figure;
scatter(x(filter),y(filter),25,'filled','MarkerFaceColor',[200 200 200]/256,'LineWidth',0.5)

hold on
errorbar(xmid, meany, stdy,'- .','CapSize',18,'MarkerSize',20,'Color','#D64D4D','LineWidth',1)

xlabel('silenced duration (min) before illumination');
ylabel('response time (min)');
xlim([time_threshold analysis_range_sil_dur])
ylim([0 analysis_range])
pbaspect([3 2 1])

%% plot response time distribution

a = response_time_final(response_time_final<=analysis_range);

response_time_fig = figure;
hold on
h = histogram(a,edges,'Normalization','probability');
h.FaceColor = mRNA_red;
mean(a)

xlim([0 analysis_range])
ylim([0 0.23])
xlabel('reactivation response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%% plot cumulative distribution

first_on_all_time_long = first_on_all_time_long(abs(last_on_all_time_long) > time_threshold);

cum_dist_fig = figure;
h = cdfplot(first_on_all_time_long);
h.LineWidth = 2;
xlabel('response time (min)')
ylabel('cumulative probability')
xlim([0 8])
ylim([0 1])

pbaspect([3 2 1])
