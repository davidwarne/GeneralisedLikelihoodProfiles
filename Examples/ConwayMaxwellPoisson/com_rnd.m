function [X] = com_rnd(n,lambda,nu)
%% COM_RND generates random variates from the  Conway-Maxwell-Poisson distribution
%  via numerical inversion of the CDF.
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
xmax_value = 99;
xmin_value = 0;
X = zeros(n,1);    
xmn = fix(xmin_value);
xmx = fix(xmax_value+1.0e-6);
y = xmn:xmx;        
m = length(y);
fy = zeros(m,1);
for j = 1:m
    yy = y(j);
    fy(j) = com_cdf(y(j),lambda,nu);
end
% numerical CDF inversion.
u = rand(n,1);	% uniform random deviate (0 : 1)
for i = 1:n    
    %Find closest "y" to urd1
    I = find(u(i) < fy);
    X(i) = y(I(1));
end    
return
