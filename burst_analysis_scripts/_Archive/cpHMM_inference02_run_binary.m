clear
close all

addpath(genpath('../lib'))
addpath(genpath('../../utilities'))

% Load data
readPath = ['../data/burst_analysis_data/'];
load([readPath 'spot_struct.mat'],'spot_struct')
load('data_ref.mat')

% set write path
dt = strrep(datestr(datetime),':','-');
writePath = [readPath dt filesep 'cpHMM_results_binary' filesep];
mkdir(writePath)

sampleSize = 3000;
nBoots = 25;
inferenceOptions.sampleSize = sampleSize;
inferenceOptions.nBoots = nBoots;
inferenceOptions.nStates = 2;
inferenceOptions.Tres = 20;
inferenceOptions.nSteps = 7;
inferenceOptions.alpha = 1.4;
inferenceOptions.eps = 1e-4;
inferenceOptions.maxWorkers = 25;
inferenceOptions.n_localEM = 25; % set num local runs
inferenceOptions.nStepsMax = 500; % set max stinferenceOptions.eps per inference

% get AP range info
ap_bounds = data_ref.ap_bounds;%[-0.02 0.02];
ap_axis = ap_bounds(1:end-1) + diff(ap_bounds)/2;

% get basic indexing info
knirps_vec = [spot_struct.meanKnirps];
ap_vec = [spot_struct.meanAP];
n_dp_vec = [spot_struct.nDP_cpHMM];
use_flags = [spot_struct.useFlag]&ap_vec>=ap_bounds(1)&ap_vec<=ap_bounds(2) & n_dp_vec >= 2*inferenceOptions.nSteps;
use_indices = find(use_flags);

% get project and experiment type vectors
embryo_id_vec = [spot_struct.masterID];

% get knirps range
knirps_bins = data_ref.knirps_bins_cpHMM;

% group data by average Knirps concentration
knirps_group_vec = discretize(knirps_vec(use_flags),knirps_bins);
knirps_group_index = unique(knirps_group_vec);

iter_i = 1;
wb = waitbar(0,'Conduction cpHMM inference...');

knirps_id_vec = [knirps_group_index(1) knirps_group_index(end)]; % compare high to low

iter_total = length(knirps_group_vec)*nBoots;
%%
for k = 1:length(knirps_id_vec)
  
    group_i = knirps_id_vec(k);
    
    % Extract subset of traces relevant to this subgroup       
    candidate_indices = use_indices(knirps_group_vec==group_i);
    candidate_embryo_vec = embryo_id_vec(candidate_indices);
    candidate_embryo_index = unique(candidate_embryo_vec);
        
    % iterate through bootstraps
    for b = 1:nBoots
        waitbar(iter_i/iter_total,wb);
        iter_i = iter_i + 1;
        
        % Initialize structures to store results
        local_struct_temp = struct;  % stores results from individual iterations
        output = struct;  % structure that will be saved to file               

        % first draw bootstrap sample of embryos     
        boot_embryos = randsample(candidate_embryo_index,length(candidate_embryo_index),true);
        [GC,GR] = groupcounts(boot_embryos');         
        emb_weight_vec = zeros(size(candidate_embryo_vec));
        for e = 1:length(GR)
            emb_weight_vec(GR(e)==candidate_embryo_vec) = GC(e);
        end
        
        % initialize cell to store inference traces
        fluo_data = {};
        kni_data = {};
        time_data = {};
        
        % take bootstrap sample             
        sample_ids = [];                                           
        ndp = 0;
        % randomly draw traces        
        while ndp < sampleSize
            tr_id = randsample(candidate_indices,1,true,emb_weight_vec);
            sample_ids = [sample_ids tr_id];
            fluo_data{end+1} = spot_struct(tr_id).fluoInterp;
            kni_data{end+1} = spot_struct(tr_id).knirps_vec_cpHMM;
            time_data{end+1} = spot_struct(tr_id).timeInterp;
            ndp = ndp + length(fluo_data{end});
        end
        
        
        % add them to data cells
        sample_particles = [spot_struct(sample_ids).masterID];
        
        % Random initialization of model parameters
        param_init = initialize_random(inferenceOptions.nStates, inferenceOptions.nSteps, fluo_data);

        % Approximate inference assuming iid data for param initialization                
        local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                            param_init.noise,inferenceOptions.nStates, inferenceOptions.nSteps, inferenceOptions.alpha, 500, 1e-4);

        noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
        v_iid = exp(local_iid_out.v_logs);  

        % create parallel pool if one does not already exist
        p = gcp('nocreate');
        if isempty(p)
            parpool(inferenceOptions.maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
        elseif p.NumWorkers > inferenceOptions.maxWorkers
            delete(gcp('nocreate')); % if pool with too many workers, delete and restart
            parpool(inferenceOptions.maxWorkers);
        end

        % conduct cpHMM inference
        parfor i_local = 1:inferenceOptions.n_localEM % Parallel Local EM 

            % Random initialization of model parameters
            param_init = initialize_random_with_priors(inferenceOptions.nStates, noise_iid, v_iid);

            % Get Intial Values
            pi0_log_init = log(param_init.pi0);
            A_log_init = log(param_init.A);
            v_init = param_init.v;                        
            noise_init = param_init.noise;

            %--------------------LocalEM Call-------------------------%
            local_out = local_em_MS2_reduced_memory(fluo_data, ...
                  v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
              inferenceOptions.alpha, inferenceOptions.nStepsMax, inferenceOptions.eps);  

            %---------------------------------------------------------%                
            % Save Results                 
            local_struct_temp(i_local).subset_id = i_local;
            local_struct_temp(i_local).logL = local_out.logL;                
            local_struct_temp(i_local).A = exp(local_out.A_log);
            local_struct_temp(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
            local_struct_temp(i_local).r = exp(local_out.v_logs).*local_out.v_signs / inferenceOptions.Tres;                                
            local_struct_temp(i_local).noise = 1/exp(local_out.lambda_log);
            local_struct_temp(i_local).pi0 = exp(local_out.pi0_log);
            local_struct_temp(i_local).total_time = local_out.n_iter;               
            local_struct_temp(i_local).soft_struct = local_out.soft_struct;               
        end

        % Record output
        [maxL, max_index] = max([local_struct_temp.logL]); % Get index of best result  

        % Save parameters from most likely local run
        output.pi0 =local_struct_temp(max_index).pi0;                        
        output.r = local_struct_temp(max_index).r(:);          
        output.noise = local_struct_temp(max_index).noise;
        output.A = local_struct_temp(max_index).A(:);
        output.A_mat = local_struct_temp(max_index).A;  
        output.max_logL = maxL;
        output.logL_results = [local_struct_temp.logL];

        % get soft-decoded structure
        output.soft_struct = local_struct_temp(max_index).soft_struct;                                                                     

        % other inference characteristics                        
        output.timeBin = 1;
%         output.apBin = inference_struct(i).ap_id;
        output.additionalBin = k;
        output.groupID = group_i;                                        

        output.truncInference = 0;
        output.iter_id = b;                        
        output.particle_ids = sample_particles;    
%         output.particle_sub_ids = particle_sub_id_data;    

        output.N = ndp;

        % save inference data used
        output.fluo_data = fluo_data;
        output.kni_mean = nanmean([kni_data{:}]);
        output.time_data = time_data;

        % Determine unique filename and sace

        % Generate filenames            
        fName_sub = ['hmm_results_group' sprintf('%03d',group_i) '_rep'];

        % generate random string
        rand_string = strrep(num2str(randsample(1:9,5,true)),' ','');
        
        % save
        out_file = [writePath fName_sub rand_string];          
        save([out_file '.mat'], 'output');           
    end        
end

delete(wb);