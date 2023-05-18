%% Example of generalised likelihood profiling using the Conway-Maxwell-Poisson 
% model as an example. This example provides a method of comparison with
% the true likelihood profile due to an accurate normalisation constant 
% approximation being available.
%
% authors: Christopher Drovandi (c.drovandi@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
%          David Warne (david.warne@qut.edu.au)
%          School of Mathematical Sciences,
%          Faculty of Science,
%          Queensland University of Technology
%
%

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

%% generate data
% true parameters 
theta_true = [4,2];
y = com_rnd(2000,4,2)
n = length(y);
rng(1337,'simdTwister');

%% Perform optimisation for observed data and numerical 
DataSim = @(theta) com_rnd(n,theta(1),theta(2)); % TODO: implement based on cmp.h util.h z.h ... seems pretty ok for the sake of consistency in code
lossfunc = @(theta,D) discrete_fisher_divergence(theta,D,@com_pdf_unnorm);

theta0 = [1,1];
lb = [1,1];
ub = [exp(3),exp(2)];
thetabnds = [lb;ub];

gridn=100;
deln = 1000;
Res = [gridn,deln];
K = 100;

%% try the generic function version
result = generalisedLikelihoodProfiles(lossfunc,y,DataSim,theta0,thetabnds,Res,K,0.05,true,true);

%% Extract output
fmin_ny        =result.fmin_ny;
theta_ny       =result.theta_ny;
theta_nyk      =result.theta_nyk;
fmin_prof      =result.fmin_prof;
theta_grid     =result.theta_grid;
delta_opt      =result.delta_opt;
delta_grid     =result.delta_grid;
lowerCI        =result.lowerCI;
upperCI        =result.upperCI;
delta_coverage =result.delta_coverage;
target_coverage=result.target_coverage;
threshold = result.threshold;

%% compure true profiles
negloglike = @(theta,y) -1.0*sum(com_logpdf(y,theta(1),theta(2)));
problem = createOptimProblem('fmincon','x0',theta0,...
        'objective', @(theta) negloglike(theta,y),'lb',lb,'ub',ub);
gs = GlobalSearch;
[theta_true_mle,theta_true_fmin] = run(gs,problem);
prof_true = zeros(size(theta_grid));
for i = 1:length(theta_grid(1,:))
    [~,prof_true(1,i)] = fmincon(@(theta) negloglike([theta_grid(1,i),theta],y),theta0(2),[],[],[],[],lb(2),ub(2));
end
for i = 1:length(theta_grid(2,:))
    [~,prof_true(2,i)] = fmincon(@(theta)negloglike([theta,theta_grid(2,i)],y),theta0(1),[],[],[],[],lb(1),ub(1));
end

%% have a look at potential bias
theta_names = {'\lambda','\nu'};
d = 2;
figure;
for j=1:d
    subplot(3,d,j);
    ksdensity(theta_nyk(:,j)','Support',[lb(j),ub(j)],'BoundaryCorrection','reflection');
    hold on
    plot(theta_ny(j),0,'rx');
    xlabel(['$\hat{\theta}_{',num2str(j),'k}$']);
    ylabel(['$p(\hat{\theta}_{',num2str(j),'k})$']); 
    legend({['density of $\hat{\theta}_{',num2str(j),',k}$'],['$\hat{\theta}_',num2str(j),'$']})
end

%% Final results plots
figure
for j = 1:d
    subplot(1,d,j)
    plot(theta_grid(j,:),exp(-delta_opt(j)*(fmin_prof(j,:)- min(fmin_prof(j,:)))),'-b');
    hold on
    plot(theta_grid(j,:),exp(-((prof_true(j,:)-min(prof_true(j,:)')))),'--r');
    hold on;
    plot([theta_true(j),theta_true(j)],[0,1],':r');
    plot([theta_ny(j),theta_ny(j)],[0,1],'-k');
    cut = - threshold;
    plot([theta_grid(j,1),theta_grid(j,end)],[exp(cut),exp(cut)],'--k')
    plot([lowerCI(j),lowerCI(j)],[min(exp(-delta_opt(j)*(fmin_prof(j,:)-fmin_ny))),max(exp(-delta_opt(j)*(fmin_prof(j,:)-fmin_ny)))],'--k');
    plot([upperCI(j),upperCI(j)],[min(exp(-delta_opt(j)*(fmin_prof(j,:)-fmin_ny))),max(exp(-delta_opt(j)*(fmin_prof(j,:)-fmin_ny)))],'--k');
    ylim([0,max(exp(-delta_opt(j)*(fmin_prof(j,:)-fmin_prof(j,:))))])
    xlim([lb(j),ub(j)]);
    xlabel(['$',theta_names{j},'$']);
    ylabel(['$\exp(-\delta DFD(\pi(\cdot | ',theta_names{j},')||\pi_n))$']);
end
