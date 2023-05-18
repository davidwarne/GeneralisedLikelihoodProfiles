function T = Stochastic_Model(LL, x0, pM, pK, pR, MC, Nm, Lbnd ,Rbnd)
%% STOCHASTIC_MDOEL Stochastic model one layer reaction diffusion advection 
%
% inputs:
% LL length of domain
% x0 release point
% pM probability of moving \in [0,1]
% pK probability of decay  \in [0,1]
% pR probability of moving right \in [0,1]
% MC number of realisations
% Nm number of raw moments
% Lbnd type of left boundary
% Rbnd type of right boundary

% authors: 
%          Matthew Simpson (matthew.simpson@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
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
T1=zeros(1,MC);

% extract boundary type user options
Lbnd_abs = strcmp(Lbnd,'absorb');
Lbnd_ref = strcmp(Lbnd,'reflect');
Rbnd_abs = strcmp(Rbnd,'absorb');
Rbnd_ref = strcmp(Rbnd,'reflect');

if Lbnd_abs && Rbnd_abs
    reflect = [];
elseif Lbnd_abs && Rbnd_ref
    reflect = LL;
elseif Lbnd_ref && Rbnd_abs
    reflect = 0;
end
nU = 10000;
Ustream = rand(nU,1);
ri = 1;

for i = 1:MC
    t=0.0;
    tau=1.0;
    x = x0;
    
%     if (x==0)
%         t=0;
%         T1(i)=t;
%     end
    
    if Lbnd_abs && Rbnd_abs
        absorbed = x <= 0 || x >= LL;
    elseif Lbnd_abs && Rbnd_ref
        absorbed = x <= 0;
    elseif Lbnd_ref && Rbnd_abs
        absorbed = x >= LL;
    end
    
    while ~absorbed
        
        R = Ustream(ri);
        S = Ustream(ri+1); ri = ri + 2;
        if ri >= nU
            ri = 1;
            Ustream = rand(nU,1);
        end
        
        if (x > 0 && x < LL)    
            if R >= pM % don't move
                x = x;
            elseif (R < pM && R < pR) % move right
                x = x + 1.0;
            elseif (R < pM && R >= pR) % move left
                x = x - 1.0;
            end
%             t = t+tau;
        end
        
        if (x==reflect)
            %R = Ustream(ri); ri = ri*(ri < nU) + 1;
            if R >= pM % don't move
                x = x;
            elseif (R < pM && R < pR) % reflected
                x = x;
            elseif (R < pM && R >= pR) 
                if reflect == 0
                    x = x + 1.0; % move right
                elseif reflect == LL
                    x = x - 1.0; % move left
                end
            end
        end
        
        if Lbnd_abs && Rbnd_abs
            absorbed = x <= 0 || x >= LL;
        elseif Lbnd_abs && Rbnd_ref
            absorbed = x <= 0;
        elseif Lbnd_ref && Rbnd_abs
            absorbed = x >= LL;
        end
        
        t = t+tau;
        
        if (S <= pK) % death
            break;
        end
    end
    
    T1(i)=t;
    
end

% output powers of exit times for empirical raw moment calculations
T = zeros(MC,Nm);
for n = 1:Nm
    T(:,n) = (T1.^n)';
end




