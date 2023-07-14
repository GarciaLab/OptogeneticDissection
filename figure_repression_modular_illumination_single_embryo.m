clc
clear
close all

addpath(genpath('./lib'))

%% create color map
%Green color map
cmap_green = [[247,252,253];[229,245,249];[204,236,230];[153,216,201];[102,194,164];[65,174,118];[35,139,69];[0,109,44];[0,68,27]];

old_stepNum = size(cmap_green,1);
new_stepNum = 256;

cmap_green_1 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,1)),linspace(0,1,new_stepNum));
cmap_green_2 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,2)),linspace(0,1,new_stepNum));
cmap_green_3 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,3)),linspace(0,1,new_stepNum));
cmap_green_new = [cmap_green_1' cmap_green_2' cmap_green_3']/256;


%Red color map
cmap_red = [[255,255,229];[255,247,188];[254,227,145];[254,196,79];[254,153,41];[236,112,20];[204,76,2];[153,52,4];[102,37,6]];

old_stepNum = size(cmap_red,1);
new_stepNum = 256;

cmap_red_1 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,1)),linspace(0,1,new_stepNum));
cmap_red_2 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,2)),linspace(0,1,new_stepNum));
cmap_red_3 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,3)),linspace(0,1,new_stepNum));

cmap_red_new = [cmap_red_1' cmap_red_2' cmap_red_3']/256;

%Color for lines
color_green = [38 143 75]/256;
mRNA_red = brighten([212 100 39]/256,.2);

%Color 
color_wt = [122 168 116]/256;
color_low = [191 213 151]/256;
color_high = [220 236 203]/256;

%% initialization
ap_lim = 0.02;
knirps_offset = 375698.13;

spot_traces_WT = [];
spot_traces_HIGH = [];
spot_traces_LOW = [];

knirps_traces_WT = [];
knirps_traces_HIGH = [];
knirps_traces_LOW = [];

%% plot mean fluorescence vs time (not aligned)

projectName = {'optokni_eve4+6_WT','optokni_eve4+6_ON_HIGH','optokni_eve4+6_ON_LOW'}; 

projectTraceNum = [4 3 1];

for i = 1:length(projectName)

    % load data
    load(['data' filesep 'main_analysis' filesep projectName{i} filesep 'spot_struct.mat'])
    
    if (projectName{i} == "optokni_eve4+6_WT")
        traceNum = projectTraceNum(1);
    else
        if (projectName{i} == "optokni_eve4+6_ON_HIGH")
            traceNum = projectTraceNum(2);
        else
            if (projectName{i} == "optokni_eve4+6_ON_LOW")
                traceNum = projectTraceNum(3);
            end
        end
    end

    time_orig_long = [];
    fluo_orig_long = [];
    frame_orig_long = [];
    knirps_orig_long = [];
    knirps_vec_long_raw = [];
    ap_vec_long = [];
    time_vec_long = [];
    post_turn_on_flags = [];
    post_turn_off_flags = [];
    
    mRNA_vec_long = [];

    for j = 1:length(spot_struct)

        if (spot_struct(j).TraceQCFlag == 1) && (spot_struct(j).setID == traceNum)

            % extract core vectors 
            fluo_vec_orig = spot_struct(j).fluo;
            time_vec_orig = spot_struct(j).time;
            frame_vec_orig = spot_struct(j).frames;
            knirps_vec_orig = spot_struct(j).rawNCProtein;
            ap_vec_orig = spot_struct(j).APPosNucleus;

            % calculate mean
            ever_on_vec(j) = any(~isnan(fluo_vec_orig));
            mean_ap_orig = nanmean(ap_vec_orig);
            mean_knirps(j) = nanmean(knirps_vec_orig);
            
            post_on_vec = zeros(size(ap_vec_orig));
            post_off_vec = zeros(size(ap_vec_orig));
            
            if ever_on_vec(j)
                start_i = find(~isnan(fluo_vec_orig),1);
                stop_i = find(~isnan(fluo_vec_orig),1,'last');
                if true%stop_i < length(fluo_vec)-10
                    post_off_vec(stop_i+1:end) = 1;
                end

                if start_i > 1
                    post_on_vec(start_i+1:end) = 1;
                end
            end
            
            if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
               time_orig_long = [time_orig_long time_vec_orig];
               frame_orig_long = [frame_orig_long frame_vec_orig];
               fluo_orig_long = [fluo_orig_long fluo_vec_orig];
               knirps_orig_long = [knirps_orig_long knirps_vec_orig];
            end

            all_zeros = fluo_vec_orig;
        
            all_zeros(isnan(all_zeros)) = 0;
            mRNA_vec_long = [mRNA_vec_long all_zeros];
            
            knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec_orig];
            ap_vec_long = [ap_vec_long ap_vec_orig];
            time_vec_long = [time_vec_long time_vec_orig];
            
            post_turn_on_flags = [post_turn_on_flags post_on_vec];
            post_turn_off_flags = [post_turn_off_flags post_off_vec];
            
        end

    end


    % calculate mean knirps and fraction on

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

        nBoots = 100;
        if ~isempty(fluo_orig_long_zero(time_filter_long))
            boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long_zero(time_filter_long));
            fluo_vec_mean(j) = nanmean(boot_samples_fluo);
            fluo_vec_ste(j) = std(boot_samples_fluo);
    
            boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_orig_long(time_filter_long));
            knirps_vec_mean(j) = nanmean(boot_samples_knirps);
            knirps_vec_ste(j) = std(boot_samples_knirps);

        else
            knirps_vec_mean(j) = NaN;
            knirps_vec_ste(j) = NaN;
        end

        frac_on(j) = nanmean(fluo_orig_long_binary(time_filter_long));
        frac_on_ste(j) = std(fluo_orig_long_binary(time_filter_long));

    end

    if (projectName{i} == "optokni_eve4+6_ON_HIGH")
        %knirps_vec_mean = calibrate_export_laser_custom(knirps_vec_mean);
        knirps_vec_mean(42:end) = calibrate_export_laser_custom(knirps_vec_mean(42:end));
    end

    knirps_vec_mean = knirps_vec_mean-knirps_offset;
    frac_on_movmean = movmean(frac_on,5);
    fluo_vec_movmean = movmean(fluo_vec_mean,10);
    fluo_vec_ste_movmean = movmean(fluo_vec_ste,10);

    if (projectName{i} == "optokni_eve4+6_WT")
        WT.time_vec = time_vec;
        WT.knirps_vec_mean = knirps_vec_mean;
        WT.knirps_vec_ste = knirps_vec_ste;
        WT.frac_on = frac_on;
        WT.frac_on_movmean = frac_on_movmean;
        WT.fluo_vec_mean = fluo_vec_mean;
        WT.fluo_vec_ste = fluo_vec_ste;
        WT.fluo_vec_movmean = fluo_vec_movmean;
        WT.fluo_vec_ste_movmean = fluo_vec_ste_movmean;
    else
        if (projectName{i} == "optokni_eve4+6_ON_HIGH")

            HIGH.time_vec = time_vec;
            HIGH.knirps_vec_mean = knirps_vec_mean;
            HIGH.knirps_vec_ste = knirps_vec_ste;
            HIGH.frac_on = frac_on;
            HIGH.frac_on_movmean = frac_on_movmean;
            HIGH.fluo_vec_mean = fluo_vec_mean;
            HIGH.fluo_vec_ste = fluo_vec_ste;
            HIGH.fluo_vec_movmean = fluo_vec_movmean;
            HIGH.fluo_vec_ste_movmean = fluo_vec_ste_movmean;
        else
            if (projectName{i} == "optokni_eve4+6_ON_LOW")
                LOW.time_vec = time_vec;
                LOW.knirps_vec_mean = knirps_vec_mean;
                LOW.knirps_vec_ste = knirps_vec_ste;
                LOW.frac_on = frac_on;
                LOW.frac_on_movmean = frac_on_movmean;
                LOW.fluo_vec_mean = fluo_vec_mean;
                LOW.fluo_vec_ste = fluo_vec_ste;
                LOW.fluo_vec_movmean = fluo_vec_movmean;
                LOW.fluo_vec_ste_movmean = fluo_vec_ste_movmean;
            end
        end
    end
    
    %% Look at long vectors

    %set parameters
    timeBins = 61;

    ap_bins = linspace(-0.2,0.2,31);
    ap_bins_plot = (ap_bins(1:end-1) + ap_bins(2:end))/2;
    time_bins = linspace(0,35*60,timeBins);
    time_bins_plot = (time_bins(1:end-1) + time_bins(2:end))/2/60;

    % calculate mean vectors
    knirps_vec_long = knirps_vec_long_raw - knirps_offset;

    ap_groups = discretize(ap_vec_long,ap_bins); 
    time_groups = discretize(time_vec_long,time_bins); 

    frac_inst_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
    frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
    eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

    frac_inst_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
    frac_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
    eve_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);

    frac_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
    frac_inst_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
    eve_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);


    frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
    frac_on_knirps_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
    eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);

    for t = 1:length(time_bins)
        for a = 1:length(ap_bins)
            time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
            if sum(time_window_filter) > 10
                still_on_flags_time = ~post_turn_off_flags(time_window_filter);
                inst_on_flags_time = mRNA_vec_long(time_window_filter)>0;
                
                frac_on_time_array_mean(t,a) = nanmean(still_on_flags_time);
                frac_on_time_array_num(t,a) = sum(~isnan(still_on_flags_time));
                frac_inst_on_time_array_mean(t,a) = mean(inst_on_flags_time);
                frac_inst_on_time_array_num(t,a) = length(inst_on_flags_time);

                eve_sample = mRNA_vec_long(time_window_filter);
                eve_sample(isnan(eve_sample)) = 0;
                eve_time_array_mean(t,a) = mean(eve_sample);
                eve_time_array_num(t,a) = sum(~isnan(eve_sample));
                eve_time_array_std(t,a) = std(eve_sample);
                eve_time_array_ste(t,a) = eve_time_array_std(t,a)/sqrt(eve_time_array_num(t,a));

                knirps_sample = knirps_vec_long(time_window_filter);
                knirps_sample(isnan(knirps_sample)) = 0;
                knirps_time_array_mean(t,a) = mean(knirps_sample);
                knirps_time_array_num(t,a) = sum(~isnan(knirps_sample));
                knirps_time_array_std(t,a) = std(knirps_sample);
                knirps_time_array_ste(t,a) = knirps_time_array_std(t,a)/sqrt(knirps_time_array_num(t,a));

            end
        end
    end
    
    knirps_time_array_mean(isnan(knirps_time_array_mean)) = 0;
    eve_time_array_mean(isnan(eve_time_array_mean)) = 0;
    frac_inst_on_time_array_mean(isnan(frac_inst_on_time_array_mean)) = 0;

    knirps_time_array_mean = movmean(knirps_time_array_mean,3,1);
    eve_time_array_mean = movmean(eve_time_array_mean,3,1);
    frac_inst_on_time_array_mean = movmean(frac_inst_on_time_array_mean,3,1);

    if i == 1
        WT.knirps_mean = knirps_time_array_mean;
        WT.eve_mean = eve_time_array_mean;
        WT.inst_on = frac_inst_on_time_array_mean;
    else
        if i == 2
            CONST.knirps_mean = knirps_time_array_mean;
            CONST.eve_mean = eve_time_array_mean;
            CONST.inst_on = frac_inst_on_time_array_mean;
        else
            LOW.knirps_mean = knirps_time_array_mean;
            LOW.eve_mean = eve_time_array_mean;
            LOW.inst_on = frac_inst_on_time_array_mean;
        end
    end

    num = 0;

    for j = 1:length(spot_struct)

        temp_trace = zeros(1,151);
        temp_knirps = NaN(1,151);

        ap_vec = spot_struct(j).APPosNucleus;
        ap_pos = mean(ap_vec);

        if (ap_pos>=-ap_lim) && (ap_pos<=ap_lim)
            num = num+1;
            spot_fluo = spot_struct(j).fluoInterp;
            knirps_protein = spot_struct(j).rawNCProteinInterp;

            frame_start = spot_struct(j).timeInterp(1)/20 + 1;
            frame_final = spot_struct(j).timeInterp(end)/20 + 1;

            if ~isnan(frame_start) && ~isnan(frame_final)
                temp_trace(frame_start:frame_final) = spot_fluo;
                temp_knirps(frame_start:frame_final) = knirps_protein;
                if i == 1
                    spot_traces_WT = [spot_traces_WT;temp_trace];
                    knirps_traces_WT = [knirps_traces_WT;temp_knirps];
                else
                    if i == 2
                        spot_traces_HIGH = [spot_traces_HIGH;temp_trace];
                        knirps_traces_HIGH = [knirps_traces_HIGH;temp_knirps];
                    else
                        spot_traces_LOW = [spot_traces_LOW;temp_trace];
                        knirps_traces_LOW = [knirps_traces_LOW;temp_knirps];
                    end
                end
            end
        end

    end
    
end


%% plot mean rate

mean_rate_comparison_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile
hold on

errorbar(LOW.time_vec,LOW.knirps_vec_mean,LOW.knirps_vec_ste,'Color','k','CapSize',0);
plot(LOW.time_vec,LOW.knirps_vec_mean,'-k','LineWidth',1)
scatter(LOW.time_vec,LOW.knirps_vec_mean,90,'d','MarkerFaceColor',color_low,'MarkerEdgeColor','k')

errorbar(HIGH.time_vec,HIGH.knirps_vec_mean,HIGH.knirps_vec_ste,'Color','k','CapSize',0);
plot(HIGH.time_vec,HIGH.knirps_vec_mean,'-k','LineWidth',1)
scatter(HIGH.time_vec,HIGH.knirps_vec_mean,90,'^','MarkerFaceColor',color_high,'MarkerEdgeColor','k')

errorbar(WT.time_vec,WT.knirps_vec_mean,WT.knirps_vec_ste,'Color','k','CapSize',0);
plot(WT.time_vec,WT.knirps_vec_mean,'-k','LineWidth',1)
scatter(WT.time_vec,WT.knirps_vec_mean,90,'MarkerFaceColor',color_wt,'MarkerEdgeColor','k')


xlim([8 30])
ylim([2.5E5 9.25E5])
xlabel(['time (min) into nc14'])
ylabel(['[Knirps] (AU)'])
pbaspect([3 1 1])

nexttile
hold on

errorbar(LOW.time_vec,LOW.fluo_vec_movmean,LOW.fluo_vec_ste_movmean,'Color','k','CapSize',0);
plot(LOW.time_vec,LOW.fluo_vec_movmean,'-k','LineWidth',1)
scatter(LOW.time_vec,LOW.fluo_vec_movmean,90,'d','MarkerFaceColor',color_low,'MarkerEdgeColor','k')

errorbar(HIGH.time_vec,HIGH.fluo_vec_movmean,HIGH.fluo_vec_ste_movmean,'Color','k','CapSize',0);
plot(HIGH.time_vec,HIGH.fluo_vec_movmean,'-k','LineWidth',1)
scatter(HIGH.time_vec,HIGH.fluo_vec_movmean,90,'^','MarkerFaceColor',color_high,'MarkerEdgeColor','k')

errorbar(WT.time_vec,WT.fluo_vec_movmean,WT.fluo_vec_ste_movmean,'Color','k','CapSize',0);
plot(WT.time_vec,WT.fluo_vec_movmean,'-k','LineWidth',1)
scatter(WT.time_vec,WT.fluo_vec_movmean,90,'MarkerFaceColor',color_wt,'MarkerEdgeColor','k')

xlim([8 30])
ylim([0 17E4])
xlabel(['time (min) into nc14'])
ylabel(['mean transcription rate (AU)'])
pbaspect([3 1 1])