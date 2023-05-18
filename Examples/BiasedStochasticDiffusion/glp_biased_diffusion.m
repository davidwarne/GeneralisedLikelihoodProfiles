%% Example of generalised likelihood profiling using the bias stochastic
% transport model as an example. This example demonstrates the applicability 
% of the approach for models where the likelihood is intractable, but a
% sensible loss function can be obtained. In this case the loss function 
% invloves an expression for the moments of the particle lifetime distribution
% 
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

validateTF = 0;

% set up experiments 
experimentrun = 1:3;
experiments.pK = [0.001,0.001,0.001];
experiments.pR = [0.48,0.50,0.52];
experiments.Nm = [2,2,2];
experiments.alpha = [0.99,0.95:-0.05:0.05,0.01]
experiments.K = [500,500,500];
experiments.deln = [1000,1000,1000];
experiments.gridn = [200,200,200];
experiments.obs_coverage_glp = zeros(2,length(experiments.alpha),3);

for id=experimentrun
%% Model parameters for synthetic "observed" data
L = 100;  % domain length
delta = 1;% step length
tau = 1;  % time step 
pM = 1.0; % probability of move event
pK = experiments.pK(id); % probability of death even
pR = experiments.pR(id); % probability of moving right given move event 
MC = 1e3; % number of stochastic simulations
x0 = round([L/8,L/4]); % point of interest
Nm = experiments.Nm(id);   % number of moments
Lbnd = 'absorb'; Rbnd = 'reflect';
rng(1337,'simdTwister');
%% Generate empirical moments
Data = bootstrapMoments(L, x0, pM, pK, pR, MC, Nm, Lbnd ,Rbnd);
fprintf('Stochastic moments and bootstap standard deviation done\n');

%% Numerical moments
% Computes numerical approximation to moments M_{n}(x) for n = 1,...,Nm by
% solving appropriate boundary value problem using a finite difference method
% Note: can also use Exact Moments using Laplace transform approach, but is
% typically much slower.
Nx = 501; % number of nodes
xn = linspace(0,L,Nx); % node locations
Mn = Numerical_Moments(L,pM,pK,pR,delta,tau,Nx,Nm,Lbnd,Rbnd); % Mn(i,n) is nth moment at x = xn(i).
fprintf('Numerical moments done\n');
Mf_numeric = @(theta,x0) Numerical_Moments_x0(L,pM,theta(1),theta(2),delta,tau,Nx,Nm,x0,Lbnd,Rbnd);

%% Perform optimisation for observed data and numerical 
DataSim = @(theta) bootstrapMoments(L,x0,pM,theta(1),theta(2),MC,Nm,Lbnd,Rbnd);
lossfunc = @(theta,D) moments_distance(theta,Mf_numeric,x0,D.Ms,D.Ms_std);
theta0 = [0.0001,0.46];
lb = [0.00001,0.3];
ub = [0.003,1];
thetabnds = [lb;ub];
Res = [experiments.gridn(id),experiments.deln(id)];
K = experiments.K(id);
alpha = experiments.alpha;

%% try the generic function version
[result] = generalisedLikelihoodProfiles(lossfunc,Data,DataSim,theta0,thetabnds,Res,K,alpha,true,true);

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
lower_bs_CI = result.lower_bs_CI;
upper_bs_CI = result.upper_bs_CI;

%% have a look at potential bias
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

figure;
d = 2
theta_names = {'d','r'};
 for j=1:d
     subplot(3,d,j);
     ksdensity(theta_nyk(:,j)','Support',[lb(j),ub(j)],'BoundaryCorrection','reflection');
     hold on
     plot(theta_ny(j),0,'rx');
     xlabel(['$\hat{p}_{',theta_names{j},',\rm{MGLE}}^{k}$']);
     ylabel(['$p(\hat{p}_{',theta_names{j},',\rm{MGLE}}^{k})$']); 
     legend({['density of $\hat{p}_{',theta_names{j},',\rm{MGLE}}^{k}$'],['$\hat{p}_{',theta_names{j},',\rm{MGLE}}$']})
 end
%% plot delta vs coverage results
for j=1:d
     subplot(3,d,d+j);
     plot(delta_grid,delta_coverage(j,:,end-1),'-');
     hold on;
     plot([min(delta_grid),max(delta_grid)],[target_coverage(end-1),target_coverage(end-1)],'--r');
     plot([delta_opt(j,end-1),delta_opt(j,end-1)],[0,1],'--k');
     xlabel(['$\delta_',num2str(j),'$']);
     ylabel(['coverage']);
     
end
legend({'coverage','target','$\delta^*$'})
%% Final results plot
theta_names = {'p_d','p_r'};
theta_true = [pK,pR];
for j = 1:d
     subplot(3,d,2*d+j)
     plot(theta_grid(j,:),exp(-delta_opt(j,end-1)*fmin_prof(j,:)),'-b');
     hold on;
     plot([theta_true(j),theta_true(j)],[0,1],'--r');
     plot([theta_ny(j),theta_ny(j)],[0,1],'-k');
     cut = -delta_opt(j,end-1)*fmin_ny - threshold(end-1);
     plot([theta_grid(j,1),theta_grid(j,end)],[exp(cut),exp(cut)],'-.k')
     plot([lowerCI(j,end-1),lowerCI(j,end-1)],[min(exp(-delta_opt(j,end-1)*fmin_prof(j,:))),max(exp(-delta_opt(j,end-1)*fmin_prof(j,:)))],'--k');
     plot([lower_bs_CI(j,end-1),lower_bs_CI(j,end-1)],[min(exp(-delta_opt(j,end-1)*fmin_prof(j,:))),max(exp(-delta_opt(j,end-1)*fmin_prof(j,:)))],':k');
     plot([upperCI(j,end-1),upperCI(j,end-1)],[min(exp(-delta_opt(j,end-1)*fmin_prof(j,:))),max(exp(-delta_opt(j,end-1)*fmin_prof(j,:)))],'--k');
     plot([upper_bs_CI(j,end-1),upper_bs_CI(j,end-1)],[min(exp(-delta_opt(j,end-1)*fmin_prof(j,:))),max(exp(-delta_opt(j,end-1)*fmin_prof(j,:)))],':k');
     
     ylim([0,max(exp(-delta_opt(j,end-1)*fmin_prof(j,:)))])
     xlim([lb(j),ub(j)]);
     xlabel(['$',theta_names{j},'$']);
     ylabel(['$\tilde{L}_{\delta^*}(',theta_names{j},' ; y)$']);
 end
legend({'profile','true parameter', 'MGLE', 'threshold' '$95\%$ CI (profile)','$95\%$ CI (quantile)'})

%% check 1: multiple data sets with true parameter
if validateTF ==1
d = 2
 theta_true = [pK,pR];
Data_reps = cell(K,1);
results_reps = cell(K,1);
obs_coverage_glp = zeros(d,1);
obs_coverage_qbs = zeros(d,1);
lowerCI_reps = zeros(d,K);
upperCI_reps = zeros(d,K);

for k=1:K
    Data_reps{k} = bootstrapMoments(L, x0, pM, pK, pR, MC, Nm, Lbnd ,Rbnd);
    results_reps{k} = generalisedLikelihoodProfiles(lossfunc,Data_reps{k},DataSim,theta0,thetabnds,Res,K,alpha,false,false);
end
for kk=1:length(experiments.alpha)

    % use calibrated confidence intervals for original dataset
    pos_opt_theta_k = zeros(d,1);
    pos_lowerCI_k = zeros(d,1);
    pos_upperCI_k = zeros(d,1);
    
    for j=1:d
        for k=1:K
            [~, pos_opt_theta_k(j)] = min(results_reps{k}.fmin_prof(j,:));
            cut = -delta_opt(j,kk)*results_reps{k}.fmin_ny - threshold(kk);
            % get grid locations of CI for base choice delta = 1
            [~,pos_lowerCI_k(j)] = min(abs(-delta_opt(j,kk)*results_reps{k}.fmin_prof(j,1:pos_opt_theta_k(j)) - cut));
            [~,pos_upperCI_theta_rev_k] = min(abs(-delta_opt(j,kk)*results_reps{k}.fmin_prof(j,end:-1:pos_opt_theta_k(j)) - cut));
            pos_upperCI_k(j) = length(results_reps{k}.fmin_prof(j,:)) - pos_upperCI_theta_rev_k +1;
            % map index to actual CI limits
            lowerCI_reps(j,k) = results_reps{k}.theta_grid(j,pos_lowerCI_k(j));
            upperCI_reps(j,k) = results_reps{k}.theta_grid(j,pos_upperCI_k(j));
        end
        
    end
    for j=1:d
        obs_coverage_glp(j,kk) = mean((theta_true(j) > lowerCI_reps(j,:)) & (theta_true(j) < upperCI_reps(j,:)),2);
    end
end


experiments.obs_coverage_glp(:,:,id) = obs_coverage_glp;
end
%% check 2: repeat for other parameter sets
end

if validateTF == 1
experiments.obs_coverage_glp
1-experiments.alpha
%%
figure;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
for i=1:3
    subplot(1,3,i);
plot(1-experiments.alpha,1-experiments.alpha,'--k','LineWidth',2)
hold on
plot(1-experiments.alpha,experiments.obs_coverage_glp(1,:,i),'-b','LineWidth',2)
plot(1-experiments.alpha,experiments.obs_coverage_glp(2,:,i),'-r','LineWidth',2)
xticks([0.0,0.25,0.5,0.75,1.0])
yticks([0.0,0.25,0.5,0.75,1.0])
xlabel('Target coverage');
ylabel('Observed coverage');
legend('optimal','$p_d$','$p_r$')
end
end
