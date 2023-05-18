function [Max0] = Exact_Moments_x0(L,pM,pK,pR,delta,tau,Nx,Nm,x0,Lbnd,Rbnd)
%% EXACT_MOMENTS_X0 Extracts the Moment for the initial position of interest x0
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

Mn = Exact_Moments(L,pM,pK,pR,delta,tau,Nm,Lbnd,Rbnd); % Ma(n) is exact moment M_{n}(x)
Mnx0 = zeros(length(x0),Nm)

for i = 1:length(x0)
    % At point of interest (x = x0)
    j = find(xn==x0(i));
    if isempty(j)
        j = find(xn-x0(i) > 0,'first'); % linearly interpolate if x0 is not an element of xn
        Mnx0(i,:) = Mn(j-1,:) + (Mn(j,:)-Mn(j-1,:))*(x0-xn(j-1))/(xn(j)-xn(j-1)); % nth moment at x = x0
    else
        Mnx0(i,:) = Mn(j,:); % nth moment at x = x0
    end
end
