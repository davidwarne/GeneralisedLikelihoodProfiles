function[pdf] = com_pdf(n, lambda, nu)
%% COM_PDF compute Conway-Maxwell-Poisson Probability Density Function.
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
summax = 100; 
termlim = 1e-6;
sum1 = 0;
logsum = -Inf;
for js = 1:summax
    logX = log(lambda^(js-1)) - nu*gammaln(js);
    M = max(logsum,logX);
    logsum = log((exp(logsum-M) + exp(logX-M)))+M;
end
pdf = exp(n*log(lambda) - (nu*gammaln(n+1) + logsum));
return
