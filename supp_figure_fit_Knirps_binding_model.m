clear
% close all


% specify number of binding sites
n_bs = 10;

% specify target function
H = -6.3;
Kd = 3.8;
kon0 = 3.2;
n_calc_points = 100;

% generate Knirps axis 
repressorCVec = linspace(0.01,10,n_calc_points);

% generate reference vurve
kon_ref = kon0 * repressorCVec.^H ./(repressorCVec.^H + Kd^H);

% generate helper functions for fit
% params = [eInteraction, eCoop, kb0, kon0]
kon_pd_fun = @(params) predict_kon_curve(repressorCVec,params(4),params(1),params(2),params(3),n_bs,n_calc_points);
loss_fun = @(params) kon_ref-kon_pd_fun(params);

params_fit = lsqnonlin(loss_fun,[0.05, 0.05, 1 3], [0, -1, 0 2],[1 0 10 5]);
% [kon_curve, n_bound_curve] = predict_kon_curve(repressorCVec,kon0,eInteraction,eCoop,kb0,n_bs,n_points);

kon_pd = kon_pd_fun(params_fit);

figure;
plot(repressorCVec,kon_ref)
hold on
plot(repressorCVec,kon_pd)
