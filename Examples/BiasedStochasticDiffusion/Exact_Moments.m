function M = Exact_Moments(L,pM,pK,pR,delta,tau,Nm,Lbnd,Rbnd)
%% EXACT_MOMENTS Computes exact moments M_{n}(x) for n = 1,..., Nm by using 
% Laplace transform approach of Ellery et al. (2013)
% Inputs
%   L: length of domain
%   pM: probability of moving
%   pK: probability of death
%   pR: probability of moving right (0 < pR < pM)
%   delta: step length
%   tau: step duration
%   Nm: number of moments
% Outputs
%   M: moment functions (vector of length Nm), where M(n) is M_{n}(x).
%
% authors: 
%          Elliot Carr (elliot.carr@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
%          David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology

pR = pR*pM;
D = delta^2*pM/(2*tau); % diffusivity
v = delta*((pM-pR)-pR)/tau; % drift/advection coefficient 
k = pK/tau; % decay rate

syms s m C x M

m(1) = (v + sqrt(v^2+4*D*(k+s)))/(2*D);
m(2) = (v - sqrt(v^2+4*D*(k+s)))/(2*D);

% Absorbing at x = 0 and reflecting at x = L

A = [1 1; m(1)*exp(m(1)*L) m(2)*exp(m(2)*L)];
b = [-1/(k+s); 0];

if strcmp(Lbnd,'absorb') && strcmp(Rbnd,'absorb')
    A = [1 1; exp(m(1)*L) exp(m(2)*L)];
    b = [-1/(k+s); -1/(k+s)];
elseif strcmp(Lbnd,'absorb') && strcmp(Rbnd,'reflect')
    A = [1 1; m(1)*exp(m(1)*L) m(2)*exp(m(2)*L)];
    b = [-1/(k+s); 0];
elseif strcmp(Lbnd,'reflect') && strcmp(Rbnd,'absorb')
    A = [m(1) m(2); exp(m(1)*L) exp(m(2)*L)];
    b = [0; -1/(k+s)];
end

a = A\b;

C = a(1)*exp(m(1)*x) + a(2)*exp(m(2)*x) + 1/(k+s);

for n = 1:Nm
    M(n) = n * (-1)^(n-1) * limit(diff(C,s,n-1),s,0);
end
