function sweepResults = sweep_loop_v3(sweepInfo,sweepResults)

% make directory to store temporary files
% tempSavePath = [savePath filesep sweepInfo.simType '_tempSweepFiles' filesep];
% mkdir(tempSavePath)

% initialize parallel poo

% initialize stuff for waitbar
WB = waitbar(0,'conducting parameter sweeps...');
D = parallel.pool.DataQueue;    
afterEach(D, @nUpdateWaitbar);

N = sweepInfo.nIterations;
p = 1;
    
% iterate through different param values
nIterations = sweepInfo.nIterations;

for sweep_step = 1:nIterations
    
    sweepInfoTemp = sweepInfo;
    
    % conduct WT simulations
    [sweepResults(sweep_step), sweepInfoTemp] = io_prediction_wrapper_cm(sweepInfoTemp,sweepResults(sweep_step));

    % conduct RA simulations
    sweepResults(sweep_step) = io_prediction_wrapper_ra_v2(sweepInfoTemp,sweepResults(sweep_step));       
    
    % update waitbar
    send(D, sweep_step);
end
delete(WB);


% helper function
function nUpdateWaitbar(~)
  waitbar(p/N, WB);
  p = p + 1;
end


end