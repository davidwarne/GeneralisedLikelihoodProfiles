function f = moments_distance(theta,Mf,x0, Ms, Ms_std)
%% MOMENTS_DISTANCE given model parameters compares empirical moments
% from observed data with analystical/numerical moments of the model at a set of
% positions of interest.
%
% Parameters:
% theta  - model parameters
% Mf     - function for computing moments
% Ms     - empirical moments of the data
% Ms_std - bootstrapped standard deviations of the moments from the data
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

    Mftheta = Mf(theta,x0);
    f = sum(sum(((Ms - Mftheta)./Ms_std).^2));
end
