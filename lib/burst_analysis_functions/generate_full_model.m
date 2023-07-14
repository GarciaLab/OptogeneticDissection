% Function wrapper for stochastic trace simulations in non-steady-state
% conditions
function sweepInfo = generate_full_model(sweepInfo,varargin)

% Build 2 state network
sweepInfo.RateMatrix = zeros(2,2);
sweepInfo.RateMatrix(2,1) = sweepInfo.kon;
sweepInfo.RateMatrix(1,2) = sweepInfo.koff;

% normalize
diag_flags = eye(size(sweepInfo.RateMatrix,1))==1;
sweepInfo.RateMatrix(diag_flags) = 0;
sweepInfo.RateMatrix(diag_flags) = -sum(sweepInfo.RateMatrix);

% add r2
sweepInfo.r_emission(2) = sweepInfo.r2;

% indicate which rate is tf-dependent (assume only one possible for now)
sweepInfo.tf_dependent_flags = false(2,2);
sweepInfo.tf_dependent_flags(2,1) = true;
