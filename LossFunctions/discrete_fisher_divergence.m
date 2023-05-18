function [sumf] = discrete_fisher_divergence(theta,y,pdf_unnorm) 
%% DISCRETE_FISHER_DIVERGENCE computes the Fisher divergence for a one
% dimensional discrete distributition with respect to an empirical
% distribution.
%
% Parameters:
% theta      - model parameters
% y          - samples from empircial distribution
% pdf_unnorm - un-normalised pdf of model
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
n = length(y);
f = zeros(n,1);
star = max(y);

for i = 1:n
    if y(i) == 0
        f(i) = (pdf_unnorm(theta, star)/pdf_unnorm(theta, y(i)))^2 ...
            - 2*pdf_unnorm(theta, y(i))/pdf_unnorm(theta, y(i)+1);
    else
        f(i) = (pdf_unnorm(theta, y(i)-1)/pdf_unnorm(theta, y(i)))^2 ...
            - 2*pdf_unnorm(theta, y(i))/pdf_unnorm(theta, y(i)+1);
    end
end
sumf = sum(f);
end


