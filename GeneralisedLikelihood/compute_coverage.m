function [coverage,lower_pjK,upper_pjK] = compute_coverage(delta, threshold, fmin_yK,fmin_pjK,pj_grid, theta_y)
%% COMPUTE_COVERAGE For a given threshold aiming for a particular coverage 
% interval and generalise likelihood parameter delta compute the empirical coverage 
% given K MLE and profile likelihoods based on simulated data at the MLE.
%
% Authors: David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
%          Christopher Drovandi (c.drovandi@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
 K = length(fmin_yK);

 lower_pjK = zeros(1,K);
 upper_pjK = zeros(1,K);

 for k = 1:K
    cut = -delta*fmin_yK(k) - threshold;

    [~, pos_opt_pjk] = min(fmin_pjK(:,k));

    [~, pos_lower_pj] = min(abs(-delta*fmin_pjK(1:pos_opt_pjk,k) - cut));
    lower_pjK(k) = pj_grid(pos_lower_pj);
    [~, pos_upper_pj_rev] = min(abs(-delta*fmin_pjK(end:-1:pos_opt_pjk,k) - cut));
    pos_upper_pj = length(fmin_pjK(:,k)) - pos_upper_pj_rev +1;
    upper_pjK(k) = pj_grid(pos_upper_pj);
 end

 coverage = mean((theta_y > lower_pjK) & (theta_y < upper_pjK));
