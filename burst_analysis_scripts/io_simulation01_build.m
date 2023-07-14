% script to call functions that build set of reference data fro parameter
% sweeps
clear
% close all
addpath(genpath('../lib'));

% set dataroot
dataRoot =['..' filesep 'data' filesep 'main_analysis' filesep];
inferenceRoot = ['..' filesep 'data' filesep 'burst_analysis_data' filesep '03-Mar-2022 17-41-41' filesep];

% specify write directory 
writeDir = ['..' filesep 'data' filesep 'burst_analysis_data' filesep];
mkdir(writeDir);

load('data_ref.mat')

% build reactivation set
projectNameRA = 'optokni_eve4+6_ON'; 
disp('Building reactivation reference set...')
tic
io_ref_ra = build_ra_ref_set(projectNameRA,dataRoot,data_ref);
toc
disp('Done.')

% build combined io dataset
projectNameCell = {'optokni_eve4+6_WT','optokni_eve4+6_ON_LOW','optokni_eve4+6_ON_HIGH'};
disp('Building combined io set...')
tic
io_ref_cm = build_combined_io_set(projectNameCell, dataRoot, data_ref);
toc
disp('Done.')

% estimate probabilistic detection threshold 
projectNameDetection = 'optokni_eve4+6_WT'; 
disp('Building detection probability reference set...')
tic
detection_stats = estimateDetectionThreshold(projectNameDetection, dataRoot, inferenceRoot);
toc

% save
disp('Saving...')
save([writeDir 'io_ref_cm.mat'],'io_ref_cm')
save([writeDir 'io_ref_ra.mat'],'io_ref_ra')
save([writeDir 'detection_stats.mat'],'detection_stats')
disp('Done.')

