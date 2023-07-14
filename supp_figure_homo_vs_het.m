clear all
close all
clc


%% parameters
knirps_offset = 0;

%% initialization

projectName = {'optoknirps_homo_control_embryo1','optoknirps_homo_control_embryo2','optoknirps_homo_control_embryo3','optoknirps_homo_control_embryo4'...
    'optoknirps_het_control_embryo1','optoknirps_het_control_embryo2','optoknirps_het_control_embryo3'};

%set parameters
timeBins = 21;

ap_bins = linspace(-0.2,0.2,31);
ap_bins_plot = (ap_bins(1:end-1) + ap_bins(2:end))/2;
time_bins = linspace(0,35,timeBins);
time_bins_plot = (time_bins(1:end-1)+time_bins(2:end))/2;

knirps_comparison = zeros(length(ap_bins)-1,length(projectName));

for k = 1:length(projectName)

    % load data
    load(['data' filesep 'knirps_homo_het' filesep projectName{k} filesep 'CompiledNuclei.mat'])

    ap_map = load(['data' filesep 'knirps_homo_het' filesep 'AP_map' filesep projectName{k} '_APmap.mat']);

    %%
    % initialize longform vectors for regression
    xcoord_vec_long = [];
    ycoord_vec_long = [];
    ap_vec_long = [];
    time_vec_long = [];
    knirps_vec_long_raw = [];


    for i = 1:length(CompiledNuclei)


        % extract core vectors 
        frame_vec = CompiledNuclei(i).Frames;
        time_vec = ElapsedTime(frame_vec)-ElapsedTime(nc14);
        knirps_vec = CompiledNuclei(i).FluoTimeTrace;
        xcoord_vec = CompiledNuclei(i).xPos;
        ycoord_vec = CompiledNuclei(i).yPos;

        %if time_vec(end) - time_vec(1) >=30*60

            % make average vectors        

            knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
            time_vec_long = [time_vec_long time_vec];
            xcoord_vec_long = [xcoord_vec_long xcoord_vec];
            ycoord_vec_long = [ycoord_vec_long ycoord_vec];
            for j = 1:length(xcoord_vec)
                ap_vec_long = [ap_vec_long ap_map.APmap(ceil(ycoord_vec(j)),ceil(xcoord_vec(j)))];
            end

        %end

    end


    %% calculate long vectors


    % calculate mean vectors
    knirps_vec_long = knirps_vec_long_raw - knirps_offset;

    ap_groups = discretize(ap_vec_long,ap_bins); 
    time_groups = discretize(time_vec_long,time_bins); 


    knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
    knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);

    for t = 1:length(time_bins)
        for a = 1:length(ap_bins)
            time_window_filter = time_groups==t & ap_groups==a;
            if sum(time_window_filter) > 10

                knirps_sample = knirps_vec_long(time_window_filter);
                %knirps_sample(isnan(knirps_sample)) = 0;
                knirps_time_array_mean(t,a) = mean(knirps_sample,'omitnan');
                knirps_time_array_num(t,a) = sum(~isnan(knirps_sample));
                knirps_time_array_std(t,a) = std(knirps_sample,'omitnan');
                knirps_time_array_ste(t,a) = knirps_time_array_std(t,a)/sqrt(knirps_time_array_num(t,a));


            end

        end
    end
    %knirps_time_array_mean(isnan(knirps_time_array_mean)) = 0;
    %knirps_time_array_mean = movmean(knirps_time_array_mean,5,1);

    
    knirps_comparison(:,k) = knirps_time_array_mean(18,:);
    

end

%% Final plot

knirps_homo = mean(knirps_comparison(:,1:4),2);
knirps_het = mean(knirps_comparison(:,5:7),2);
    
hold on
plot(ap_bins_plot,knirps_comparison(:,1:4),'-o','LineWidth',1.5,'Color',[0.75 0.75 0.75])
plot(ap_bins_plot,knirps_comparison(:,5:7),'-*','LineWidth',1.5,'Color',[0.75 0.75 0.75])
plot(ap_bins_plot,knirps_homo,'LineWidth',3);
plot(ap_bins_plot,knirps_het,'LineWidth',3);
xlim([-0.17 0.2])
ylim([1.75E5 16E5])
xlabel('AP position')
ylabel('fluorescence (AU)')
legend('','','','','','','','Homo', 'Het')

pbaspect([3 2 1])