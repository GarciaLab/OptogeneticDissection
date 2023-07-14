% Script to compile inference results 
clear
close all

% Load data

dateString = '21-Jul-2022 10-12-27';

dataRoot = '..\data\burst_analysis_data\';
dataPath = [dataRoot dateString filesep];
load([dataRoot 'inference_data.mat'],'inference_data')
load([dataPath 'inference_info.mat'],'inference_struct')

% set path to raw cpHMM output
inferencePath = [dataPath 'cpHMM_results*' filesep];

% get list of inference results
inf_files = dir([inferencePath 'hmm_results*']);

% read raw inference files into memory
inference_results = struct;
wb = waitbar(0,'loading inference results...');
for i = 1:length(inf_files)
    waitbar(i/length(inf_files),wb);
    load([inf_files(i).folder filesep inf_files(i).name])
    f_mean = nanmean([output.fluo_data{:}]);
    t_mean = nanmean([output.time_data{:}]);
    fnames = fieldnames(output);
    for f = 1:length(fnames)
        inference_results(i).(fnames{f}) = output.(fnames{f});
    end
    inference_results(i).f_mean = f_mean;
    inference_results(i).t_mean = t_mean;
end  
delete(wb)


z_index = [inference_struct.z_id];
particle_id_vec_long = inference_data.particle_id_vec;
group_id_vec = [inference_results.groupID];
knirps_id_vec = [inference_results.additionalBin];
gp_kni_array = unique([group_id_vec' knirps_id_vec'],'rows');
diff_vec = [1 diff(gp_kni_array(:,2))'];
inf_ids = cumsum(1*(diff_vec<0))+1;

knirps_axis_vec = [inference_results.kni_mean];

%% generate full inference id vec
inf_id_vec = NaN(size(group_id_vec));
for i = 1:length(inf_id_vec)
    ind = find(group_id_vec(i)==gp_kni_array(:,1) & knirps_id_vec(i)==gp_kni_array(:,2));
    inf_id_vec(i) = inf_ids(ind);
end
z_flag_vec = z_index(inf_id_vec);

% generate result vectors and arrays
A_array = cat(3,inference_results.A_mat);
R_array = NaN(size(A_array));
use_flags = false(size(knirps_axis_vec));

for i = 1:size(A_array,3)
   A_temp = A_array(:,:,i);
   warning('') % Clear last warning message
   [R_temp, exit_flag] = logm(A_temp);
   [warnMsg, warnId] = lastwarn;
   R_array(:,:,i) = R_temp/20;
   use_flags(i) = all(isreal(R_temp(:))) && sum(R_temp(:)>=0)==2 && isempty(warnMsg);
end  

r_array = cat(2,inference_results.r);

logL_vec = [inference_results.max_logL];

noise_vec = [inference_results.noise];

pi0_array = vertcat(inference_results.pi0);

pon_vec = reshape(A_array(2,1,:),[],1);
poff_vec = reshape(A_array(1,2,:),[],1);

kon_vec = reshape(R_array(2,1,:),[],1);
koff_vec = reshape(R_array(1,2,:),[],1);

group_index = unique(group_id_vec);
knirps_mean_vec = NaN(size(group_index));
knirps_ste_vec = NaN(size(group_index));
zeros_vec  = NaN(size(group_index));

kon_mean_vec = NaN(size(group_index));
kon_ste_vec = NaN(size(group_index));
kon_upper_vec = NaN(size(group_index));
kon_lower_vec = NaN(size(group_index));

koff_mean_vec = NaN(size(group_index));
koff_ste_vec = NaN(size(group_index));
koff_upper_vec = NaN(size(group_index));
koff_lower_vec = NaN(size(group_index));

dur_mean_vec = NaN(size(group_index));
dur_ste_vec = NaN(size(group_index));

pon_mean_vec = NaN(size(group_index));
pon_upper_vec = NaN(size(group_index));
pon_lower_vec = NaN(size(group_index));

poff_mean_vec = NaN(size(group_index));
poff_upper_vec = NaN(size(group_index));
poff_lower_vec = NaN(size(group_index));

r2_mean_vec = NaN(size(group_index));
r2_ste_vec = NaN(size(group_index));
r2_upper_vec = NaN(size(group_index));
r2_lower_vec = NaN(size(group_index));

r1_mean_vec = NaN(size(group_index));
r1_ste_vec = NaN(size(group_index));
r1_upper_vec = NaN(size(group_index));
r1_lower_vec = NaN(size(group_index));

pi0_mean_array = NaN(size(group_index,1),length(inference_results(1).pi0));
noise_mean_vec = NaN(size(group_index));
noise_ste_vec = NaN(size(group_index));
noise_upper_vec = NaN(size(group_index));
noise_lower_vec = NaN(size(group_index));

for g = 1:length(group_index)
  
    grp_indices = find(group_id_vec==group_index(g)&use_flags);
    if length(grp_indices) > 5
        % check for outliers in the logL values
        [~, outlier_flags] = rmoutliers(logL_vec(grp_indices));

        knirps_mean_vec(g) = nanmean(knirps_axis_vec(grp_indices(~outlier_flags)));
        knirps_ste_vec(g) = nanstd(knirps_axis_vec(grp_indices(~outlier_flags)));

        zeros_vec(g) = unique(z_flag_vec(grp_indices(~outlier_flags)));

        kon_mean_vec(g) = nanmean(kon_vec(grp_indices(~outlier_flags)));
        kon_ste_vec(g) = nanstd(kon_vec(grp_indices(~outlier_flags)));
        kon_upper_vec(g) = prctile(kon_vec(grp_indices(~outlier_flags)),75);
        kon_lower_vec(g) = prctile(kon_vec(grp_indices(~outlier_flags)),25);

        noise_mean_vec(g) = nanmean(noise_vec(grp_indices(~outlier_flags)));
        noise_ste_vec(g) = nanstd(noise_vec(grp_indices(~outlier_flags)));
        noise_upper_vec(g) = prctile(noise_vec(grp_indices(~outlier_flags)),75);
        noise_lower_vec(g) = prctile(noise_vec(grp_indices(~outlier_flags)),25);

        pi0_mean_array(g,:) = mean(pi0_array(grp_indices(~outlier_flags),:),1);

        koff_mean_vec(g) = nanmean(koff_vec(grp_indices(~outlier_flags)));
        koff_ste_vec(g) = nanstd(koff_vec(grp_indices(~outlier_flags)));
        koff_upper_vec(g) = prctile(koff_vec(grp_indices(~outlier_flags)),75);%nanstd(koff_vec(grp_indices(~outlier_flags)));
        koff_lower_vec(g) = prctile(koff_vec(grp_indices(~outlier_flags)),25);

        dur_mean_vec(g) = nanmean(1./koff_vec(grp_indices(~outlier_flags)));
        dur_ste_vec(g) = nanstd(1./koff_vec(grp_indices(~outlier_flags)));

        pon_mean_vec(g) = nanmean(pon_vec(grp_indices(~outlier_flags)));
        pon_upper_vec(g) = prctile(pon_vec(grp_indices(~outlier_flags)),75);%nanstd(pon_vec(grp_indices(~outlier_flags)));
        pon_lower_vec(g) = prctile(pon_vec(grp_indices(~outlier_flags)),25);

        poff_mean_vec(g) = nanmean(poff_vec(grp_indices(~outlier_flags)));
        poff_upper_vec(g) = prctile(poff_vec(grp_indices(~outlier_flags)),75);%nanstd(poff_vec(grp_indices(~outlier_flags)));
        poff_lower_vec(g) = prctile(poff_vec(grp_indices(~outlier_flags)),25);

        r2_mean_vec(g) = nanmean(r_array(2,grp_indices(~outlier_flags)));
        r2_ste_vec(g) = nanstd(r_array(2,grp_indices(~outlier_flags)));
        r2_upper_vec(g) = prctile(r_array(2,grp_indices(~outlier_flags)),75); %r_ste_vec(g) = nanstd(r_array(2,grp_indices(~outlier_flags)));
        r2_lower_vec(g) = prctile(r_array(2,grp_indices(~outlier_flags)),25);

        r1_mean_vec(g) = nanmean(r_array(1,grp_indices(~outlier_flags)));
        r1_ste_vec(g) = nanstd(r_array(1,grp_indices(~outlier_flags)));
        r1_upper_vec(g) = prctile(r_array(1,grp_indices(~outlier_flags)),75); %r_ste_vec(g) = nanstd(r_array(2,grp_indices(~outlier_flags)));
        r1_lower_vec(g) = prctile(r_array(1,grp_indices(~outlier_flags)),25);
    end
end    

% Compile results into grouped vectors and save as matlab structure

%% generate indexing vectors
inf_id_index = [1];% 2 3];
zeros_flags = [0];% 0 0];
% exp_type_cell = {'ON LOW', 'Wildtype', 'ON HIGH'};
inference_summary = struct;
iter = 1;
for ind = inf_id_index
  
    % obtain filtering vector
    plot_filter = inf_ids==ind;
    % NL: excluding one obvious outlier
    if ind == 3
        last_i = find(plot_filter,1,'last')
        plot_filter(last_i) = false;
    end
    % save metadata
%     inference_summary(iter).exp_type = exp_type_cell{iter};
    inference_summary(iter).zero_flag = zeros_flags(iter);
    inference_summary(iter).inf_indices = find(plot_filter);
    inference_summary(iter).inf_sub = find(ismember(group_id_vec,group_index(plot_filter)));
    
    % save results
    inference_summary(iter).knirps_ste = knirps_ste_vec(plot_filter); 
    inference_summary(iter).knirps_mean = knirps_mean_vec(plot_filter);     
        
    inference_summary(iter).kon_upper = kon_upper_vec(plot_filter)*60; 
    inference_summary(iter).kon_lower = kon_lower_vec(plot_filter)*60; 
    inference_summary(iter).kon_ste = kon_ste_vec(plot_filter)*60; 
    inference_summary(iter).kon_mean = kon_mean_vec(plot_filter)*60; 

    inference_summary(iter).pi0_mean = pi0_mean_array(plot_filter,:); 
    
    inference_summary(iter).koff_upper = koff_upper_vec(plot_filter)*60; 
    inference_summary(iter).koff_lower = koff_lower_vec(plot_filter)*60; 
    inference_summary(iter).koff_ste = koff_ste_vec(plot_filter)*60; 
    inference_summary(iter).koff_mean = koff_mean_vec(plot_filter)*60; 
    
    inference_summary(iter).dur_ste = dur_ste_vec(plot_filter)*60; 
    inference_summary(iter).dur_mean = dur_mean_vec(plot_filter)*60; 
    
    inference_summary(iter).pon_upper = pon_upper_vec(plot_filter); 
    inference_summary(iter).pon_lower = pon_lower_vec(plot_filter); 
    inference_summary(iter).pon_mean = pon_mean_vec(plot_filter); 
    
    inference_summary(iter).poff_upper = poff_upper_vec(plot_filter); 
    inference_summary(iter).poff_lower = poff_lower_vec(plot_filter); 
    inference_summary(iter).poff_mean = poff_mean_vec(plot_filter); 
    
    inference_summary(iter).r2_upper = r2_upper_vec(plot_filter); 
    inference_summary(iter).r2_lower = r2_lower_vec(plot_filter); 
    inference_summary(iter).r2_ste = r2_ste_vec(plot_filter); 
    inference_summary(iter).r2_mean = r2_mean_vec(plot_filter); 
    
    inference_summary(iter).r1_upper = r1_upper_vec(plot_filter); 
    inference_summary(iter).r1_lower = r1_lower_vec(plot_filter); 
    inference_summary(iter).r1_ste = r1_ste_vec(plot_filter); 
    inference_summary(iter).r1_mean = r1_mean_vec(plot_filter); 
    
    inference_summary(iter).noise_upper = sqrt(noise_upper_vec(plot_filter)); 
    inference_summary(iter).noise_lower = sqrt(noise_lower_vec(plot_filter)); 
    inference_summary(iter).noise_ste = sqrt(noise_ste_vec(plot_filter)); 
    inference_summary(iter).noise_mean = sqrt(noise_mean_vec(plot_filter)); 
    
    iter = iter + 1;
end  


save([dataPath filesep 'inference_summary.mat'],'inference_summary')
