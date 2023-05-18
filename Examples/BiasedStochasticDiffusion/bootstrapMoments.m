function [Retval] = bootstrapMoments(L, x0, pM, pK, pR, MC, Nm, Lbnd ,Rbnd)
%% BOOTSTRAPMOMENTS Non-parametric bootstrap for the variability of the
% Monte Carlo estimator of the particle lifetime raw moments. 
% authors: 
%          David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
%          Christopher Drovandi (c.drovandi@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%


Ms = zeros(length(x0),Nm);
Ms_std = zeros(length(x0),Nm);
B = 1000; % for bootstrap
for i = 1:length(x0)
    T = Stochastic_Model(L,x0(i),pM,pK,pR,MC,Nm,Lbnd,Rbnd); % T(:,n) is 
    for n = 1:Nm
        Ms(i,n) = mean(T(:,n)); % nth moment at x = xn(i)
    end
    % bootstrap standard deviations of moments
    Msb = zeros(B,Nm);
    for b = 1:B
        r = randsample(1:MC,MC,'true');
        Tb = T(r,:);
        for n=1:Nm
            Msb(b,:) = mean(Tb(:,n));
        end
    end
    Ms_std(i,:) = std(Msb);
end

Retval = struct();
Retval.Ms = Ms;
Retval.Ms_std = Ms_std;

end

