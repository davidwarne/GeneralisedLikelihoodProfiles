function p = pdf_unnorm(theta,x)
%% COM_PDF_UNNORM compute Conway-Maxwell-Poisson Unnormalised Probability 
% Density Function.
%
% auhtor: 
%          David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
p = (theta(1)^x)*factorial(x)^(-theta(2));
end
