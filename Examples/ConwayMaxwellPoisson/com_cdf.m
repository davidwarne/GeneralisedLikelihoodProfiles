function [cdf] = com_cdf(n, lambda, nu)
%% COM_CDF compute Conway-Maxwell-Poisson Cumulative Distribution Function/
%
%   See "Conjugate Analysis of the Conway-Maxwell-Poisson Distribution", 
%   J. Kadane et al., Carnegie Mellon et al., 6/27/20031.
%
% Based on the implementation of James Huntley (2023).
%     Generation of Random Variates 
%     (https://www.mathworks.com/matlabcentral/fileexchange/35008-generation-of-random-variates), 
%      MATLAB Central File Exchange. Retrieved 14 April, 2023.
%
% modified by: 
%          David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
sum1 = com_pdf(0, lambda, nu);
for jn = 1:n
    sum1 = sum1 + com_pdf(jn, lambda, nu);
end
cdf = sum1;
return

