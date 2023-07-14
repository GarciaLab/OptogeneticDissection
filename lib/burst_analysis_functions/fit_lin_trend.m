function [lin_trend, lin_95, lin_05, trend_array, param_95, param_05, param_array] = ...
                                                        fit_lin_trend(x_vec_out,x_vec_in,...
                                                        x_vec_ste,y_vec,y_vec_ste,weight_vec,nBoots)

    trend_array = NaN(length(x_vec_out),nBoots); 
    param_array = NaN(2,nBoots);
    for n = 1:nBoots
%         s_ids = randsample(1:length(y_vec),length(y_vec),true);    
        s_ids = 1:length(y_vec);
        
        % draw samples assuming gaussian errors
        x_vec_in = normrnd(x_vec_in(s_ids),x_vec_ste(s_ids));
        y_vec_in = normrnd(y_vec(s_ids),y_vec_ste(s_ids)); 
        
        % define error-weighted function
        mdl = fitlm(x_vec_in', y_vec_in','Weights',weight_vec);    
        trend_array(:,n) = predict(mdl,x_vec_out');  
        param_array(:,n) = mdl.Coefficients.Estimate;
    end

    % mdl_full_dur = fitlm([knirps_vec_long' (knirps_vec_long').^2],y,'Weights',weight_vec);
    mdl_full = fitlm(x_vec_in,y_vec,'Weights',weight_vec);

    % dur_trend_mean = predict(mdl_full_dur,[x' (x').^2]);
    lin_trend = predict(mdl_full,x_vec_out');
    lin_95 = prctile(trend_array,95,2)'; % NL: might want to revisit this
    lin_05 = prctile(trend_array,5,2)';
    
    % identify param percentiles
%     [~, norm_ind] = min(abs(x_vec_out-3));
    
    % 95th 
%     trend_95 = prctile(trend_array(norm_ind,:),95);
%     [~,i_95] = min(abs(trend_array(norm_ind,:)-trend_95));
    param_95 = prctile(param_array,95,2);
    
    % 5th 
%     trend_05 = prctile(trend_array(norm_ind,:),5);
%     [~,i_05] = min(abs(trend_array(norm_ind,:)-trend_05));
    param_05 = prctile(param_array,5,2);       