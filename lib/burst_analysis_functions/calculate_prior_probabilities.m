function   logL_prior = calculate_prior_probabilities(prior_mean_vec, prior_std_vec, new_params_prop)

    prior_probs = normpdf(new_params_prop(1:7), prior_mean_vec, prior_std_vec) .* sqrt(2*pi*prior_std_vec.^2);   
    logL_prior = sum(log(prior_probs));%log(p_HC) + log(p_KD) + log(p_kon) + log(p_koff) + log(p_r2);