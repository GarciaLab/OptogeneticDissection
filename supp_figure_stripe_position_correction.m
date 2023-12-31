
% script to correct stripe position based on the curvature

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))

%% initialization
% Specify the dataset to analyze

% We need to specify the prefix of the project we wish to analyze


Prefix = 'ap_correction'; APflip = 1;
avrInit = 30;
avrFinal = 35;
analysisInit = 0;
analysisFinal = 38;

%% basic settings

timeFinal = 0:1/3:45;

% initialize binning parameters
binNum = 26;
edges = linspace(-0.2,0.2,binNum);

%% Part 1: Add working path and load the data

% Now, load all the data

% Now we'll load the CompiledParticles dataset
% This generates a string that gives the path to CompiledParticles
LoadParticlePath = ['data' filesep Prefix filesep '/CompiledParticles.mat']; 
LoadNucleiPath = ['data' filesep Prefix filesep Prefix '_lin.mat'];
LoadComNucleiPath = ['data' filesep Prefix filesep 'CompiledNuclei.mat']; 
LoadAPDetPath = ['data' filesep Prefix filesep 'APDetection.mat'];

data.nuclei = load(LoadComNucleiPath);
data.sch = load(LoadNucleiPath);
data.particle = load(LoadParticlePath);
APDetData = load(LoadAPDetPath);

ElapsedTime = data.nuclei.ElapsedTime-data.nuclei.ElapsedTime(data.nuclei.nc14);

% find frames to analyze
initFrame = find(ElapsedTime<=analysisInit,1,'last');
finalFrame = find(ElapsedTime>analysisFinal,1);

% total number of frames
totalFrame = finalFrame-initFrame + 1;

% find frames to average
avrInitFrame = find(ElapsedTime<=avrInit,1,'last');
avrFinalFrame = find(ElapsedTime>avrFinal,1);

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
processed_data(1).SpotFluo = [];


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
            processed_data(frameTemp).SpotFluo = [processed_data(frameTemp).SpotFluo, 0];
            processed_data(frameTemp).NuclearFluo = [processed_data(frameTemp).NuclearFluo fluo];  
        catch
            processed_data(frameTemp).xcoord = x_coord;
            processed_data(frameTemp).ycoord = y_coord;
            processed_data(frameTemp).schnitznum = schTemp;
            processed_data(frameTemp).SpotFluo = 0;
            processed_data(frameTemp).NuclearFluo = fluo;
        end
    end
end

CompiledParticles = data.particle.CompiledParticles;

% Assign particle info to all the nuclei
for i = 1:size(CompiledParticles{1,1},2)
    schnitz_num = CompiledParticles{1,1}(i).schnitz;
    for j = 1:size(CompiledParticles{1,1}(i).Frame,2)
        frame = CompiledParticles{1,1}(i).Frame(j);
        if frame<=finalFrame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            processed_data(frame).SpotFluo(num) = CompiledParticles{1,1}(i).Fluo(j);
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
        CompiledNucleiData(i).spotFluo(frameTemp) = processed_data(j).SpotFluo(i);
        
    end
end
% Calculate average coordinate for each nuclei

for i = 1:nucleiNum
    CompiledNucleiData(i).xcoordMean = mean(CompiledNucleiData(i).xcoord);
    CompiledNucleiData(i).ycoordMean = mean(CompiledNucleiData(i).ycoord);
end

%% Part 4: Fit the protein average for a moving time window

x_coord = zeros(length(schPass),1);
y_coord = zeros(length(schPass),1);
fluo = zeros(length(schPass),1);

% Compile all the nuclei in each frame and assign basic info
for i = 1:length(schPass)

    sch_num = schPass(i);
    index_s = find(schnitzcells(sch_num).frames == avrInitFrame);
    index_f = find(schnitzcells(sch_num).frames == avrFinalFrame);

    x_coord(i) = mean(schnitzcells(sch_num).cenx(index_s:index_f));
    y_coord(i) = mean(schnitzcells(sch_num).ceny(index_s:index_f));
    fluo(i) = mean(max(schnitzcells(sch_num).Fluo(index_s:index_f,:)));

end

%%

sf1 = fit([x_coord,y_coord],fluo,'loess'); % good, local quadratic regression

NuclearFitFig = figure;
plot(sf1,[x_coord,y_coord],fluo);
xlabel('x position (pixels)')
ylabel('y position (pixels)')
zlabel('Knirps concentration (AU)')
daspect([1 1 3000])

%% Part 5: Calibrate AP position
% Correction for the curvature
% Try fitting the center line to quadratic function

% Calculate coarse grained image
xscale = 5;
yscale = 5;

xlimit = 0:xscale:768;
ylimit = 0:yscale:450;

f  = @(x,y)sf1(x,y)-1E5;

[X Y] = meshgrid(xlimit,ylimit);
F = f(X,Y);

maxPeak = max(max(F));
% Extract stripe region
ImMask = F>maxPeak*0.8;
ImFluo = F.*ImMask;

ImFluo_rot = imrotate(ImFluo,-90);
factor = length(xlimit);
k_temp = find(ImFluo_rot);
x_temp = floor(k_temp/factor)*xscale+1;
y_temp = mod(k_temp,factor)*yscale;
w = ImFluo_rot(k_temp);

%define model
modelFun2 = @(b,x) b(1)*x.^2+b(2)*x+b(3);
%modelFun3 = @(b,x) b(1)*x.^3+b(2)*x.^2+b(3)*x+b(4);

%initial condition
start2 = [0 0 0];
%start3 = [0 0 0 0];

%fit without weight
% %nlm3 = fitnlm(X(:),Y(:),modelFun2,start2);
%fit with weight
%wnlm3 = fitnlm(x_temp,y_temp,modelFun3,start3,'Weight',w);
wnlm2 = fitnlm(x_temp,y_temp,modelFun2,start2,'Weight',w);

fig = figure;
x_temp = [0 450];
y_temp = [0 768];
%imagesc(x_temp,y_temp,ImFluo_rot);
imagesc(y_temp,x_temp,ImFluo_rot');
xx = linspace(-200,650)';

%plot stripe
line(predict(wnlm2,xx),xx,'color','b','LineWidth',3)
ylim([0 450])
xlim([0 768])
xlabel('x position (pixels)')
ylabel('y position (pixels)')
axis equal
% Calculate distance to the central curve

% initialize a meshgrid and central curve
xscale = 1;
yscale = 1;

xlimit = 1:xscale:768;
ylimit = 1:yscale:450;
[X Y] = meshgrid(xlimit,ylimit);

point_all = [X(:) Y(:)];
curve_all = [predict(wnlm2,xx) xx];

% calculate distance of each point to the central curve
[xy,distance,t] = distance2curve(curve_all,point_all);
point_predict = predict(wnlm2,point_all(:,2));
flag = point_predict>point_all(:,1);
distance(flag) = -distance(flag);

APmap = flip(reshape(distance,450,768),1);

if APflip
    APmap = -APmap;
end

% Correction for embryo length

APLength = sqrt(sum((APDetData.coordAZoom-APDetData.coordPZoom).^2));

APmap = APmap/APLength;

fig = figure;
imagesc(y_temp,450-x_temp,APmap*100)
xx = linspace(-200,650)';
%plot stripe
line(predict(wnlm2,xx),xx,'color','b','LineWidth',3)

xlabel('x position (pixels)')
ylabel('y position (pixels)')
ylim([0 450])
xlim([0 768])

axis equal
%lcolorbar('distance from center (% embryo length)')
colorbar
colormap(jet)
caxis([-20 20])
