function M = Numerical_Moments(L,pM,pK,pR,delta,tau,Nx,Nm,Lbnd,Rbnd)
%% NUMERICAL_MOMENTS Computes numerical approximation to moments M_{n}(x) for 
% n = 1,...,Nm by solving appropriate boundary value problem using a finite 
% difference method
% Inputs
%   L: length of domain
%   pM: probability of moving
%   pK: probability of death
%   pR: probability of moving right (0 < pR < pM)
%   delta: step length
%   tau: step duration
%   Nm: number of moments
% Outputs
%   M: moments (matrix of size Nx by Nm), where M(i,n) is M_{n}(x(i)),
%   where x(i) = (i-1)*L/(Nx-1)
%
% authors: 
%          Elliot Carr (elliot.carre@qut.edu.au)
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

% Finite difference discretisation
h = L/(Nx-1); % node spacing
%A = zeros(Nx,Nx); 
b = -ones(Nx,1);

% direct construnction of sparse matrix
e = ones(Nx,1);
A = spdiags([(D/h^2 + v/(2*h))*e, (-2*D/h^2 - k)*e, (D/h^2 - v/(2*h))*e],-1:1,Nx,Nx);

% left and right boundary
if strcmp(Lbnd,'absorb')
    A(1,1) = 1; b(1) = 0; % absorbing boundary
elseif strcmp(Lbnd,'reflect')
    A(1,1) = -2*D/h^2 - k; A(1,2) = 2*D/h^2;
end

if strcmp(Rbnd,'absorb')
    A(Nx,Nx) = 1; b(Nx) = 0; % absorbing boundary
elseif strcmp(Rbnd,'reflect')
    A(Nx,Nx-1) = 2*D/h^2; A(Nx,Nx) = -2*D/h^2 - k; % reflecting boundary
end
b = sparse(b);

%  for i = 2:Nx-1
%      A(i,i-1) = D/h^2 + v/(2*h);
%      A(i,i) = -2*D/h^2 - k;
%      A(i,i+1) = D/h^2 - v/(2*h);
%  end
% 
% A = sparse(A);
% 
% Moments
M = zeros(Nx,Nm);
M(:,1) = A\b; % moment 1
for n = 2:Nm % moments 2,...,Nm
    b = -n*M(:,n-1);
    if strcmp(Lbnd,'absorb')
        b(1) = 0;
    end
    if strcmp(Rbnd,'absorb')
        b(Nx) = 0;
    end
    M(:,n) = A\b;
end
