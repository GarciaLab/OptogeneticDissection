function ra_cdf_pd_mean = calculate_ra_cdf(off_flags, before_flags, fluo_array, index_vec, max_time, sweepInfo, thresh_val)

    % ra_flags 
    ra_indices = find(off_flags&before_flags);    
    reactivation_time_vec = NaN(1,size(fluo_array,2));
    
    % set all leading and trailing zeros to NaN
    for i = ra_indices        
        fluo_vec = fluo_array(:,i);        
        
        % find last "OFF" frame
        temp = fluo_vec;
        temp(index_vec<0) = thresh_val;        
        last_i = find(temp>thresh_val,1); 
        
        if isempty(last_i)
            reactivation_time_vec(i) = max_time+100;
        else
            switch_time = index_vec(last_i);
            % estimate reactivation time                    
            reactivation_time_vec(i) = switch_time*sweepInfo.deltaT;  
        end                                                                                             
    end
            
    % construct empirical cdf for ractivation
    ra_time_vec = sweepInfo.reactivation_time;

    % conduct bootstrapping                
    ra_time_vec_nn = reactivation_time_vec(~isnan(reactivation_time_vec));
    ra_count_raw = (0:length(ra_time_vec_nn))/length(ra_time_vec_nn); 
    
    % sort   
    [ra_times_sorted,~] = sort(ra_time_vec_nn);
    
    % generate dummy time vector so that matlab won't throw a fit
    bs_time_vec = sort([0 ra_times_sorted+rand(size(ra_times_sorted))*1e-6]);      
    if length(bs_time_vec) > 1
        try
            if max(bs_time_vec) < ra_time_vec(end)
                bs_time_vec(end+1) = ra_time_vec(end);
                ra_count_raw(end+1) = ra_count_raw(end);
            end
            ra_cdf_pd_mean = interp1(bs_time_vec,ra_count_raw,ra_time_vec,'previous');    
        catch
            ra_cdf_pd_mean = zeros(size(ra_time_vec));    
        end
    else
        ra_cdf_pd_mean = zeros(size(ra_time_vec));    
    end   