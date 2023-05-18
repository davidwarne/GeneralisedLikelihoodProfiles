function [result] = generalisedLikelihoodProfiles(lossfunc,Data,sim,theta0,thetabnds,Res,K,alpha,gsTF,calibrateTF)
% GENERALISEDLIKELIHOODPROFILES Implements the estimation of the
% Maximum Generalised Likelihood for a simulation-based inference problem.
% Univariate Generalised Likelihood profiles are constructed and rescaled
% for a target freqentist coverage level. 
%
% Inputs:
%  lossfunc    - loss function with signature @(theta,Data) where theta is the 
%                parameter vector and Data is the observed data. must be strictly 
%                positive function.
%  Data        - object containing observed data.
%  sim         - stochastic simulation function that generated realisations of the 
%                assumed data generating process. Signature is @(theta) where theta
%                are the model parameters.
%  theta0      - initial parameter estimates for optimisation routines.
%  thetabnds   - Bound of on parameter space with thetabnds(1,:) being lower bounds`
%                and thetabnds(2,:) being upper bounds. 
%  Res         - vector of grid resolutions, Res(1) is the profile grid resolution
%                and Res(2) is the resolution for optimising delta for calibration.
%  K           - number of parametric bootstrap samples to perfomr
%  alpha       - vector of traget significance levels
%  gsTF        - toggle global search for MGLE (true = global search)
%  calibrateTF - toggle calibration of tuning parameter targing coverage 1 - alpha
%                if false then return raw unscaled profiles with delta = 1.
%
% Output:
% result - a struct containing at least the following elements:
%
%    result.fmin_ny          - loss function at the MGLE.
%    result.theta_ny         - the MGLE.
%    result.fmin_prof        - unscaled loss function profile.
%    result.theta_grid       - grids for each profile.
%
% if calibrateTF = true then these addition elements are included.
%
%    result.theta_nyk        - K bootstrap MGLEs.
%    result.delta_opt        - optimal tuning parameter for each profile at each
%                              target coverage level.
%    result.delta_grid       - grid for with delta was optimised over;
%    result.lowerCI          - profile based lower confidence bound for each 
%                              parameter at each target coverage level.
%    result.upperCI          - profile based upper confidence bound for each 
%                              parameter at each target coverage level.
%    result.delta_coverage   - empirical coverage obtained during calibration.
%    result.target_coverage  - target coverage list (just 1 - alpha).
%    result.threshold        - list of thresholds for each confidence level;
%    result.lower_bs_CI      - quantile bootstrap based lower confidence bound 
%                              for each parameter at each target coverage level.
%    result.upper_bs_CI      - quantile bootstrap based upper confidence bound 
%                              for each parameter at each target coverage level.
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


gridn = Res(1);
deln = Res(2);
d = length(theta0);
lb = thetabnds(1,:);
ub = thetabnds(2,:);

%% Minimise loss function to obtain MGLE

% set up minimisation problem
problem = createOptimProblem('fmincon','x0',theta0,...
        'objective', @(theta)lossfunc(theta,Data),'lb',lb,'ub',ub);
% optimise    
if gsTF == true
    gs = GlobalSearch; % run global optimisation sweep, needed for some
                       % complex loss surfaces
    [theta_ny,fmin_ny] = run(gs,problem);
else
    [theta_ny,fmin_ny] = fmincon(problem);
end

%% Generate unscaled generalised likelihood profiles for each parameter
theta_grid = zeros(d,gridn);
fmin_prof = zeros(size(theta_grid));
for j=1:d
    theta_grid(j,:) = linspace(lb(j),ub(j),gridn);
    for i=1:gridn
        [~,fmin_prof(j,i)] = fmincon(@(theta)lossfunc([theta(1:(j-1)),theta_grid(j,i),theta(j:(d-1))],Data),...
            theta_ny([1:(j-1),(j+1):d]),[],[],[],[],lb([1:(j-1),(j+1):d]),ub([1:(j-1),(j+1):d]));
    end
end


%% Generate K calibration simulations
if calibrateTF == true
    fmin_profK = zeros(d,gridn,K);
    fmin_nyk = zeros(1,K);
    theta_nyk = zeros(K,d);

    parfor k = 1:K
        % simulated dataset given MGLE
        Data_k = sim(theta_ny);
        % MGLEs and profiles
        options = optimoptions('fmincon','Display','notify-detailed');
        [theta_nyk_tmp,fmin_nyk(k)] = fmincon(@(theta)lossfunc(theta,Data_k),...
            theta0,[],[],[],[],lb,ub,[],options);
        for j=1:d
            for i=1:gridn
                [~,fmin_profK(j,i,k)] = fmincon(@(theta)lossfunc([theta(1:(j-1)),theta_grid(j,i),theta(j:(d-1))],Data_k),...
                theta_nyk_tmp([1:(j-1),(j+1):d]),[],[],[],[],lb([1:(j-1),(j+1):d]),ub([1:(j-1),(j+1):d]));
            end
        end
        theta_nyk(k,:) = theta_nyk_tmp;
    end
 
    %% optimise scaling for target coverage

    % get the appropriate threshold
    ad = length(alpha);
    target_coverage = 1- alpha;
    threshold = chi2inv(target_coverage,1)/2;

    delta_grid = linspace(0.0,10,deln); 
    delta_opt = ones(d,ad);
    delta_coverage = zeros(d, deln,ad);
    theta_coverage = zeros(d,ad);
    
    for k=1:ad
        for i=1:deln
            for j=1:d
                [delta_coverage(j,i,k)] = compute_coverage(delta_grid(i),threshold(k),fmin_nyk, ...
                    squeeze(fmin_profK(j,:,:)),theta_grid(j,:),theta_ny(j));
                if abs(theta_coverage(j,k) - target_coverage(k)) > abs(delta_coverage(j,i,k) - target_coverage(k))
                    delta_opt(j,k) = delta_grid(i);
                    theta_coverage(j,k) = delta_coverage(j,i,k);
                end
            end
        end
    end

    % estimate Confidence intervals for original dataset
    pos_opt_theta = zeros(d,1);
    pos_lowerCI = zeros(d,1);
    pos_upperCI = zeros(d,1);
    lowerCI = zeros(d,ad);
    upperCI = zeros(d,ad);
    for k=1:ad
        for j=1:d
            [~, pos_opt_theta(j)] = min(fmin_prof(j,:));
            cut = -delta_opt(j,k)*fmin_ny - threshold(k);
            % get grid locations of CI 
            [~,pos_lowerCI(j)] = min(abs(-delta_opt(j,k)*fmin_prof(j,1:pos_opt_theta(j)) - cut));
            [~,pos_upperCI_theta_rev] = min(abs(-delta_opt(j,k)*fmin_prof(j,end:-1:pos_opt_theta(j)) - cut));
            pos_upperCI(j) = length(fmin_prof(j,:)) - pos_upperCI_theta_rev +1;
            % map index to actual CI limits
            lowerCI(j,k) = theta_grid(j,pos_lowerCI(j));
            upperCI(j,k) = theta_grid(j,pos_upperCI(j));
        end
    end
    % because we already have the bootstrap samples, compute typical bootstrap
    % CIs
    lower_bs_CI = zeros(d,ad);
    upper_bs_CI = zeros(d,ad);
    for k=1:ad
        for j=1:d
            lower_bs_CI(j,k) = 2*theta_ny(j) - quantile(theta_nyk(:,j),1-alpha(k)/2);
            upper_bs_CI(j,k) = 2*theta_ny(j) - quantile(theta_nyk(:,j),alpha(k)/2);
        end
    end
    
    result = struct();
    result.fmin_ny          = fmin_ny;
    result.theta_ny         = theta_ny;
    result.fmin_prof        = fmin_prof;
    result.theta_grid       = theta_grid;
    result.theta_nyk        = theta_nyk;
    result.delta_opt        = delta_opt;
    result.delta_grid       = delta_grid;
    result.lowerCI          = lowerCI;
    result.upperCI          = upperCI;
    result.delta_coverage   = delta_coverage;
    result.target_coverage  = target_coverage;
    result.threshold        = threshold;
    result.lower_bs_CI      = lower_bs_CI;
    result.upper_bs_CI      = upper_bs_CI;
else
    % output MGLE and GL profile for delta = 1
    result = struct();
    result.fmin_ny          = fmin_ny;
    result.theta_ny         = theta_ny;
    result.fmin_prof        = fmin_prof;
    result.theta_grid       = theta_grid;
end
