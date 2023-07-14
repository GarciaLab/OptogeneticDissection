
clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))

%% initialize parameters

on_frame = 43;
off_frame = 84;

analysisInit = 0;
analysisFinal = 17.5;

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

eYFP_background = 375698.13;

fontSize = 10;

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%Color 
color_wt = [122 168 116]/256;
color_low = [191 213 151]/256;
color_high = [220 236 203]/256;

%% Part 1: Add working path and load the data
%% 
% Now, load all the data

% Now we'll load the CompiledParticles dataset
% This generates a string that gives the path to CompiledParticles
LoadNucleiPath = ['data' filesep 'optogenetics_parameter' filesep 'optoknirps_eve4_6_import_export' filesep 'optoknirps_eve4_6_import_export_lin.mat'];
LoadComNucleiPath = ['data' filesep 'optogenetics_parameter' filesep 'optoknirps_eve4_6_import_export' filesep 'CompiledNuclei.mat']; 
LoadWTPath = ['data' filesep 'optogenetics_parameter' filesep 'optoknirps_eve4_6_wildtype' filesep 'data.mat']; 

data.nuclei = load(LoadComNucleiPath);
data.sch = load(LoadNucleiPath);
data.WT = load(LoadWTPath);

ElapsedTime = data.nuclei.ElapsedTime-data.nuclei.ElapsedTime(data.nuclei.nc14);

% find frames to analyze
initFrame = find(ElapsedTime<=analysisInit,1,'last');
finalFrame = find(ElapsedTime>analysisFinal,1);

% total number of frames
totalFrame = finalFrame-initFrame + 1;

totalTime = ElapsedTime(initFrame:finalFrame);
%% Part 2: Quality check on nuclei tracking

% Quality check and save the nuclei that passed the test
% Requires the nuclei to have continuous trace until lateral movement
X_pass = [];
Y_pass = [];

X_fail = [];
Y_fail = [];

schPass = [];

schnitzcells = data.sch.schnitzcells;

for i = 1:size(schnitzcells,2)
    index_s = find(schnitzcells(i).frames == initFrame);
    index_f = find(schnitzcells(i).frames == finalFrame);
    fluo_temp = max(schnitzcells(i).Fluo(index_s:index_f,:),[],2);
    result = sum(isnan(fluo_temp));
    % quality check
    if ~isempty(index_s) && ~isempty(index_f) && (index_f-index_s == finalFrame-initFrame)  ...
            && (result == 0)
        schPass = [schPass, i];
        X_pass = [X_pass, schnitzcells(i).cenx(index_s)];
        Y_pass = [Y_pass, schnitzcells(i).ceny(index_s)];
    else
        X_fail = [X_fail, schnitzcells(i).cenx(index_s)];
        Y_fail = [Y_fail, schnitzcells(i).ceny(index_s)];
    end
end


%% Part 3: Compile the schnitzcells together
% compile based on each frame

nucleiNum = length(schPass);

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];

processed_data(1).schnitznum = [];
processed_data(1).NuclearFluo = [];


% Compile all the nuclei in each frame and assign basic info
for i = 1:nucleiNum
    schTemp = schPass(i);
    index_s = find(schnitzcells(schTemp).frames == initFrame);
    index_f = find(schnitzcells(schTemp).frames == finalFrame);
    
    for j = index_s:index_f
        frameTemp = schnitzcells(schTemp).frames(j);
        x_coord = schnitzcells(schTemp).cenx(j);
        y_coord = schnitzcells(schTemp).ceny(j);
        fluo = max(schnitzcells(schTemp).Fluo(j,:));
        try
            processed_data(frameTemp).xcoord = [processed_data(frameTemp).xcoord, x_coord];
            processed_data(frameTemp).ycoord = [processed_data(frameTemp).ycoord, y_coord];
            processed_data(frameTemp).schnitznum = [processed_data(frameTemp).schnitznum, schTemp];
            processed_data(frameTemp).NuclearFluo = [processed_data(frameTemp).NuclearFluo fluo];  
        catch
            processed_data(frameTemp).xcoord = x_coord;
            processed_data(frameTemp).ycoord = y_coord;
            processed_data(frameTemp).schnitznum = schTemp;
            processed_data(frameTemp).NuclearFluo = fluo;
        end
    end
end

% Compile based on each nuclei

for i = 1:nucleiNum
    for j = initFrame:finalFrame
        
        % Transform the actual frame count to analyzed frame count
        frameTemp = j-initFrame+1;
        
        CompiledNucleiData(i).xcoord(frameTemp) = processed_data(j).xcoord(i);
        CompiledNucleiData(i).ycoord(frameTemp) = processed_data(j).ycoord(i);
        
        CompiledNucleiData(i).schnitzNum(frameTemp) = processed_data(j).schnitznum(i);
        CompiledNucleiData(i).nuclearFluo(frameTemp) = processed_data(j).NuclearFluo(i);
        
    end
end
% Calculate average coordinate for each nuclei

for i = 1:nucleiNum
    CompiledNucleiData(i).xcoordMean = mean(CompiledNucleiData(i).xcoord);
    CompiledNucleiData(i).ycoordMean = mean(CompiledNucleiData(i).ycoord);
end


%% Correct AP position

ap_map = load(['data' filesep 'optogenetics_parameter' filesep 'optoknirps_eve4_6_import_export' filesep 'optoknirps_eve4_6_import_export_APmap.mat']);

for i = 1:nucleiNum
    CompiledNucleiData(i).APPos = ap_map.APmap(ceil(CompiledNucleiData(i).ycoordMean),ceil(CompiledNucleiData(i).xcoordMean));
end

%% Calculate mean

nucleiList = [];

for i = 1:nucleiNum
    
    if (CompiledNucleiData(i).APPos>-ap_lim) & (CompiledNucleiData(i).APPos<ap_lim)
       nucleiList = [nucleiList i]; 
    end
    
end

knirps_vec_mean = zeros(totalFrame,1);
knirps_vec_ste = zeros(totalFrame,1);

for i = 1:totalFrame
    
    fluo_temp = [];
    
    for j = 1:length(nucleiList)  
        fluo_temp = [fluo_temp CompiledNucleiData(nucleiList(j)).nuclearFluo(i)];
    end
    
    nBoots = 100;
    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_temp);
    knirps_vec_mean(i) = nanmean(boot_samples_fluo);
    knirps_vec_ste(i) = std(boot_samples_fluo,'omitnan');

    if i==80
        temp = std(boot_samples_fluo,'omitnan');
    end

    %knirps_vec_mean(i) = nanmean(fluo_temp);
    %knirps_vec_ste(i) = nanstd(fluo_temp);
    
end

%% Plot result

time_interp = -2:1/10:12;

time_vec = ElapsedTime(initFrame:finalFrame)-ElapsedTime(on_frame);

knirps_vec_mean(on_frame:off_frame) = calibrate_export_laser_custom(knirps_vec_mean(on_frame:off_frame));
knirps_vec_mean = knirps_vec_mean-eYFP_background;

knirps_vec_mean_interp = interp1(time_vec,knirps_vec_mean,time_interp);

knirps_vec_ste_interp = interp1(time_vec,knirps_vec_ste,time_interp);

time_shifted = data.WT.time_plot-15.5;

start_frame_opto = find(time_interp==0);
start_frame_WT = find(time_shifted==0);

knirps_WT_rescaled = data.WT.knirps_vec_mean/data.WT.knirps_vec_mean(start_frame_WT);
knirps_WT_rescaled_ste = data.WT.knirps_vec_ste/data.WT.knirps_vec_mean(start_frame_WT);

knirps_opto_rescaled = knirps_vec_mean_interp/knirps_vec_mean_interp(start_frame_opto);
knirps_opto_ste_rescaled = knirps_vec_ste_interp/knirps_vec_mean_interp(start_frame_opto);

fig = figure;

hold on

errorbar(time_interp,knirps_vec_mean_interp,knirps_vec_ste_interp,'Color','k','CapSize',0);
plot(time_interp,knirps_vec_mean_interp,'-k','LineWidth',1)
scatter(time_interp,knirps_vec_mean_interp,60,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
xlim([-2 10])
%ylim([0.4 1.4])
ylim([3.5E5 11E5])
%ylim([3.2E5 11E5])

pbaspect([1.6 1 1])