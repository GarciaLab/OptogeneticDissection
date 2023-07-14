clear
close all

addpath(genpath('../lib'))

% Load data
readPath = '..\data\burst_analysis_data\';

% load spot data
load([readPath 'spot_struct.mat'],'spot_struct')

% load cpHMM results
dateString = '21-Jul-2022 10-12-27';
cpHMMPath = [readPath dateString filesep];
load([cpHMMPath filesep 'inference_summary.mat'],'inference_summary')

% extract relenvat parameters for Viterbi fits
tres = spot_struct(1).tresInterp;
inf_ind = 1;
v = [inference_summary.r1_mean(inf_ind) inference_summary.r2_mean(inf_ind)]*tres;
K = 2;
w = 7;
kappa = 1.4;
noise = inference_summary.noise_mean(inf_ind);
pi0_log = log(inference_summary.pi0_mean(inf_ind,:)/sum(inference_summary.pi0_mean(inf_ind,:)));
Q = [0 inference_summary.koff_mean(inf_ind) ; inference_summary.kon_mean(inf_ind) 0]/60;
Q(eye(K)==1) = -sum(Q,1);
A = expm(Q*tres);
A_log = log(A);

%%
minDP = w;
NumWorkers = 24;
pool = gcp('nocreate');
if isempty(pool)
  parpool(NumWorkers);  
elseif  pool.NumWorkers ~= NumWorkers     
  delete(pool)
  parpool(NumWorkers);  
end  
parfor i = 1:length(spot_struct)
    if sum(~isnan(spot_struct(i).fluoInterp)) >= minDP
        
      f_vec_i = spot_struct(i).fluoInterp;
      t_vec_i = spot_struct(i).timeInterp;
      t_vec_long = t_vec_i(1):tres:ceil(spot_struct(i).time(end)/tres)*tres;
      f_vec_long = zeros(size(t_vec_long));
      f_vec_long(ismember(t_vec_long,t_vec_i)) = f_vec_i;

      viterbi_out = viterbi(f_vec_long, v, noise, pi0_log, ...
                                           A_log, K, w, kappa);

      % add to spot struct
      spot_struct(i).z_viterbi_interp = viterbi_out.z_viterbi;
      spot_struct(i).t_viterbi_interp = t_vec_long;
      % adjust it to have same frame spacing as raw data
      zi = viterbi_out.z_viterbi;
      z_raw = interp1(t_vec_long,double(zi),spot_struct(i).time,'nearest','extrap');
      spot_struct(i).z_viterbi_raw = z_raw;
      
    else
      spot_struct(i).z_viterbi_raw = NaN;
      spot_struct(i).t_viterbi_interp = NaN;
      spot_struct(i).z_viterbi_interp = NaN;
    end
end

save([readPath 'spot_struct.mat'],'spot_struct')