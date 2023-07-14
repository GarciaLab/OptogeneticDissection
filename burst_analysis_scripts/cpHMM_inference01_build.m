clear
close all

addpath(genpath('../lib'))
% addpath(genpath('../../utilities'))

% Load data
readPath = '..\data\main_analysis\';
savePath = '..\data\burst_analysis_data\';

% indicate projects to use
projectNameCell = {'optokni_eve4+6_ON_LOW','optokni_eve4+6_WT','optokni_eve4+6_ON_HIGH'};
load('data_ref');

% specify correction parameters for blue light
knirps_offset = data_ref.knirps_offset; % NL: would be good to import this from elsewhere
min_dp = data_ref.min_dp;
% time_bounds = data_ref.time_bounds; % nuclei must have been around for full extent of this interval 
% ap_bounds = data_ref.ap_bounds;
knirps_cal_slope = data_ref.knirps_cal_slope;
knirps_cal_intercept = data_ref.knirps_cal_intercept;
knirps_norm_factor = data_ref.knirps_norm_factor;
% knirps_offset = 375000 / 1e5;
% knirps_cal_slope = 1.243;
% knirps_cal_intercept = 1.079e5 / 1e5; %NL: dividing everything through by 1e5 for simplicity

% load data
spot_struct_full = [];
for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};    
    load([readPath projectName filesep 'spot_struct.mat'],'spot_struct')
    for s = 1:length(spot_struct)
        spot_struct(s).projectName = projectName;
        spot_struct(s).projectID = p;
    end
    spot_struct_full = [spot_struct_full spot_struct];
end

spot_struct = spot_struct_full;

% generate a new set of unique identifiers
set_vec = [spot_struct.setID];
project_vec = [spot_struct.projectID];
project_set_array = unique([[spot_struct.projectID]' [spot_struct.setID]'],'rows');

for i = 1:size(project_set_array,1)
    projectID = project_set_array(i,1);
    setID = project_set_array(i,2);
    temp_ids = find(project_vec==projectID&set_vec==setID);
    for t = temp_ids
       spot_struct(t).masterID = i;
       particleID = spot_struct(t).particleID;
       particleIDNew = particleID - floor(particleID) + i;
       spot_struct(t).particleID = particleIDNew;
       spot_struct(t).particleIDOrig = particleID;
    end
end    

% identify and correct for effects of blue light laser
minDP = 20;
spot_struct_temp = spot_struct;
master_id_vec = [spot_struct.masterID];
master_id_index = unique(master_id_vec);
time_index_interp = 0:spot_struct(1).tresInterp:(50*60);

% NL: obtained these frames via manual inspection of protein trends
blue_light_frame_vec = [NaN(1,8) 41 33 36];
% blue_light_frame_vec = [52 NaN 52 73 46 NaN 50 36 38 miNaN 42 43 NaN 67 NaN 45 43];% OG

% initialize array to store mean protein trend
mean_protein_array = NaN(length(time_index_interp),length(master_id_index));

for m = 1:length(master_id_index)
    master_ids = find(master_id_vec==master_id_index(m));
    time_index = unique(round([spot_struct(master_ids).time],0));
    projectName = spot_struct(master_ids(1)).projectName;
                      
    % find changepoint  
    shift_frame = blue_light_frame_vec(m);
    if ~isnan(shift_frame)
        shift_time = time_index(shift_frame);               
    else
        shift_time = Inf;
    end
    % check that corrections look reasonable
    protein_array_temp = NaN(length(time_index),length(master_ids));        
    for i = master_ids
        t_vec = round(spot_struct(i).time,0);
        t_vec_fluo = spot_struct(i).timeInterp;
        start_i = find(time_index_interp<=t_vec(1)&time_index_interp<=t_vec_fluo(1),1,'last');
        stop_i = find(time_index_interp>=t_vec(end)&time_index_interp>=t_vec_fluo(end),1);
        t_vec_interp = time_index_interp(start_i:stop_i);           

        % extract basic vectors
        pt_vec = spot_struct(i).rawNCProtein/knirps_norm_factor;
        pt_vec_orig = spot_struct(i).rawNCProteinInterp/knirps_norm_factor;
        fluo_vec = spot_struct(i).fluoInterp;        
        ap_vec = spot_struct(i).APPosNucleus;
        % update time vector            
%         spot_struct(i).timeInterpOrig = spot_struct(i).timeInterp;
        spot_struct(i).timeNew = t_vec_interp;
        
        if length(t_vec) >= minDP && ~any(isnan(pt_vec)) && sum(~isnan(t_vec_interp)) >= minDP/2
            spot_struct(i).useFlag = true;
            spot_struct(i).fluoFlag = false;
            spot_struct(i).nDP_cpHMM = length(t_vec_fluo);
            raw_adjusted = pt_vec;
            raw_adjusted_orig = pt_vec_orig;
            pert_ind = find(t_vec>=shift_time,1);
            pert_ind_orig = find(t_vec_fluo>=shift_time,1);
%             pert_ind = max([pert_ind,pert_ind-1])
            if ~isempty(pert_ind) 
                raw_adjusted(pert_ind:end) = ...
                (raw_adjusted(pert_ind:end)-knirps_cal_intercept)/knirps_cal_slope;
                raw_adjusted_orig(pert_ind_orig:end) = ...
                (raw_adjusted_orig(pert_ind_orig:end)-knirps_cal_intercept)/knirps_cal_slope;
            end
            raw_adjusted = raw_adjusted - knirps_offset;
            raw_adjusted_orig = raw_adjusted_orig - knirps_offset;
            pt_interp = interp1(t_vec,raw_adjusted,t_vec_interp,'linear','extrap');
            ap_interp = interp1(t_vec,ap_vec,t_vec_interp,'linear','extrap');
            spot_struct(i).rawNCProteinNew = pt_interp;
            spot_struct(i).knirps_vec_cpHMM = raw_adjusted_orig; 
            spot_struct(i).meanKnirps = mean(raw_adjusted_orig);
            spot_struct(i).fluoNew = NaN(size(pt_interp));
            spot_struct(i).apPosNucleusNew = ap_interp;
            spot_struct(i).meanAP = mean(ap_interp(ismember(t_vec_interp,t_vec_fluo)));
%             
%             if spot_struct(i).meanKnirps < 0 
%                 error('wtf')
%             end
            if sum(~isnan(t_vec_fluo)) >= minDP/2
                spot_struct(i).fluoFlag = true;
                fluo_new = zeros(size(t_vec_interp));
                fluo_new(ismember(t_vec_interp,t_vec_fluo)) = fluo_vec;           
                spot_struct(i).fluoNew = fluo_new;
            end
        else
            spot_struct(i).useFlag = false;
            spot_struct(i).fluoFlag = false;
            spot_struct(i).rawNCProteinNew = NaN(size(t_vec_interp));
            spot_struct(i).apPosNucleusNew = NaN(size(t_vec_interp));
            spot_struct(i).fluoNew = NaN(size(t_vec_interp));
            spot_struct(i).knirps_vec_cpHMM = NaN(size(t_vec_fluo)); 
            spot_struct(i).meanKnirps = NaN;
            spot_struct(i).meanAP = NaN;
            spot_struct(i).nDP_cpHMM = NaN;
        end    
    end        

    % check that corrections look reasonable
    temp_array = NaN(length(time_index_interp),length(master_ids));    
    iter = 1;
    for i = master_ids
        t_vec = spot_struct(i).timeNew;
        if spot_struct(i).useFlag
            temp_array(ismember(time_index_interp,t_vec),iter) = spot_struct(i).rawNCProteinNew;
        end
        iter = iter + 1;
    end
    mean_protein_array(:,m) = nanmean(temp_array,2);    
end    


% Chop traces up into 15 minute pieces to try and extend the dynamic range

shift_inc = round(1*60 / spot_struct(1).tresInterp);
window_size = round(15*60 / spot_struct(1).tresInterp);
ap_bounds = [-0.15 0.15];

inf_ind_vec = find([spot_struct.fluoFlag]);
inference_data = struct;

% initialize key values vectors
inference_data.fluo_vec_cell = {};
inference_data.knirps_vec_cell = {};
inference_data.time_vec_cell = {};
inference_data.ap_vec_cell = {};
inference_data.mean_fluo_vec = [];
inference_data.mean_knirps_vec = [];
inference_data.mean_time_vec = [];
inference_data.mean_ap_vec = [];

% initialize key metadata vectors
inference_data.project_id_vec = [];
inference_data.particle_id_vec = [];
inference_data.rep_id_vec = [];
inference_data.zero_flag_vec = [];

iter = 1;
wb = waitbar(0,'building inference set');
for i = inf_ind_vec
   waitbar(iter/length(inf_ind_vec),wb);
   
   % meta data
   particleID = spot_struct(i).particleID;
   projectID = spot_struct(i).projectID;
   
   % values
   fluo_vec = spot_struct(i).fluoNew;
   first_i = find(fluo_vec>0,1)+6; % start 6 steps in because inference assumes fully loaded gene to start
   fluo_vec = fluo_vec(first_i:end);
   time_vec = spot_struct(i).timeNew(first_i:end);   
   kni_vec = spot_struct(i).rawNCProteinNew(first_i:end);
   ap_vec = spot_struct(i).apPosNucleusNew(first_i:end);
   
   % calculat how many distinct shifts we have
   n_shifts = 1 + floor((length(time_vec)-window_size)/shift_inc);
   last_obs = find(fluo_vec>0,1,'last');
   for n = 1:n_shifts
      start_i = shift_inc*(n-1)+1;
      last_i = start_i + window_size-1;
      ap_temp = ap_vec(start_i:last_i);
      ap_mean = nanmean(ap_temp);
      
      if ap_mean >= ap_bounds(1) && ap_mean <= ap_bounds(2) && ~isempty(last_obs)
          % record values
          inference_data.fluo_vec_cell{iter} = fluo_vec(start_i:last_i);
          inference_data.knirps_vec_cell{iter} = kni_vec(start_i:last_i);
          inference_data.time_vec_cell{iter} = time_vec(start_i:last_i);
          inference_data.ap_vec_cell{iter} = ap_vec(start_i:last_i);
          inference_data.mean_fluo_vec(iter) = nanmean(fluo_vec(start_i:last_i));
          inference_data.mean_knirps_vec(iter) = nanmean(kni_vec(start_i:last_i));
          inference_data.mean_time_vec(iter) = nanmean(time_vec(start_i:last_i));
          inference_data.mean_ap_vec(iter) = ap_mean;
          % meta-data
          inference_data.particle_id_vec(iter) = particleID;
          inference_data.project_id_vec(iter) = projectID;
          inference_data.rep_id_vec(iter) = n;
          inference_data.zero_flag_vec(iter) = last_i>last_obs; % signals that trace contains trailing zeros
          inference_data.zero_count_vec(iter) = max([0 last_i-last_obs]);
          % increment
          iter = iter + 1;
      end     
   end
end  
delete(wb)

save([savePath 'inference_data.mat'],'inference_data')
save([savePath 'spot_struct.mat'],'spot_struct')
