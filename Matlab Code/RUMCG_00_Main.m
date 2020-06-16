%% Code Description: Main File
% This code creates normalized budgets and computes expenditure share on a
% aggregate commodity bundles using UK budget data from 1975-1999.
% Next, it computes matrices X and A for each budget; X is a matrix of
% valid patches and A is a matrix of pseudo-agents that are have 
% preferences that is consistent with GARP.
% Last, we bootstrap a test-statistic that determines whether or not
% individuals in the UK sample have preferences consistent with RUM+GARP.
% Under the null, households satisfy RUM+GARP.

% Original Code by:
%	Yuichi Kitamura
%	Joerg Stoye
%	Ben Friedrich       (RA)
%   Duk Gyoo Kim        (RA)
%	Matthew Thirkettle  (RA)
%   Yang Zhang          (RA)

% Reference paper: 
%	Non-parametric analysis of random utility models (RUM): testing.  
%	By Yuichi Kitamura and Joerg Stoye
%	Econometrica, 2018, pp. 1883-1909

% Adapted for Column Generation by:
%   Bart Smeulders

% Reference Paper:
%   Nonparametric Analysis of Random Utility Models: Computational Tools for Statistical Testing
%   By Bart Smeulders, Laurens Cherchye and Bram de Rock

% USER INSTRUCTIONS:
% CPLEX INSTRUCTIONS:
% (1) Request CPLEX Academic licencese through the IBM Academic Initiative        
% (2) ....

% PARAMETER INSTRUCTIONS
%   (1) budget_length is the length of sub-budgets.  For example, if
%   budget_length=8, then we compute test statistics for the time periods
%   1975-1982, 1976-1983, ect.
%
%   (2) number_classes is the number of aggregate commodity bundles.
%   For example, number_classes = 3 corresponds to the commodity bundles
%   {food, services, durables}. 
%   Currently number_classes in {3,4,5} are valid inputs.  Alter
%   RUM_21_budgets.m to change aggregation methods.
%   
%   WARNING:
%   Running time increases exponentially with budget length. Especially for
%   4 and 5 goods, 15+ periods may take significant amounts of time.
%
%   (4) polynomial_degree is the polynomial degree used in the series 
%   estimation.  
%
%   (5) estimator in {0,1,2} computes the test statistic based on 
%   kernel estimator (estimator=0), series estimator (estimator = 1)
%   or endogenous estimator (estimator =2).
%   
%   (6) genAX indicators to compute matrix of pseudo-agents  A,X and comp
%   If A,X not saved (Input folder), then program computes A,X.  
%   WARNING: If one changes the price sequence or budget set, then do not
%   use an old version of A,X, as they will be incorrect.
%  
%   (7) indJ indicator to compute the test statistic.
%
%   In practice, we should always set genAX=1 and indJ=1 ... this is mainly
%   used as a debugging device, so that we do not have to recompute AX when
%   checking to see if the test statistic works.
%
%   (8) cores pre-specifies the number of cores to use.  This is to prevent
%   overuse of server resources.
%
%   (9) bs_reps is the number of bootstrap repititions used to
%   calculate test statistic.  
%
%   Other parameters:
%   (a) start_year, end_year are both determined by the dataset.  E.g. the
%   current dataset starts at 1975 and ends in 1999.
%
%   (b) seed is specified so that a new seed is used for each combination
%   of (budget_length,number_classes) to ensure independence).

%% Clear Memory and start timer
clear
tic  

warning('off','MATLAB:nargchk:deprecated')
warning('off','MATLAB:nearlySingularMatrix')

%% User-specified parameters
budget_length    = 7;       % Number of budgets
number_classes   = 3;       % Number of aggregate commodity bundles
polynomial_degree= 3;       % Series degree estimator
estimator        = 2;       % 2 = Endogenous Series, Only endogenous series possible in this package.
genAX            = 1;       % Generate pseudo-agent matrix A and patches X. Must be re-calculated (Set to 0) if tau_set_size changes
tau_set_size     = 1000;    % Number of rational choice types used in the tightening procedure
indJ             = 1;       % Run bootstrap
cores            = 4;       % Number of cores for parallel computing 
bs_reps          = 1000;    % Number of bootstrap repitions
tau_ind          = 1;       % Set equal to 0/1 to set tau equal to 0/not 0.

%% Global parameters
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores tau_val


%% Flags
% Choose estimator.  Set:
% flag_estimator = 0 for Kernel Smoothing
% flag_estimator = 1 for Basis Smoothing
% flag_estimator = 2 for Endogenous Basis Smoothing
flag_estimator  = estimator;
tau_val         = tau_ind;

% flag_genAX = 1 iff we want to generate matrix A and X and there does not
% exist a pre-generated AX file.
flag_genAX = genAX;
nameAX= strcat('Input/RUM_',num2str(budget_length),'LengthBudget_',num2str(number_classes),'.mat');
if exist(nameAX) ~= 2 & genAX == 0
    flag_genAX = 1;
end

%% Parameters
% Length of sub-budgets to be called (e.g. budget_l = 8 implies budgets 
% span 75-82, 76-83, ect.).  When there are 3 classes use budget_l <9,
% otherwise it takes a long time to run.
budget_l    = budget_length;     

% The number of commodity bundles is called classes e.g. 3 is food, 
% services, durables).  Model currently allows for classes = 3, 4, or 5.
% This can be further disaggregated by modifying RUM_21_budgets.m
classes     = number_classes;      

% Time parameters:
start_year  = 75;                        % First year for budgets
end_year    = 99;                        % Last year for budgets
periods     = end_year - start_year + 1; % Total number of periods
budget_n    = periods - (budget_l-1);    % Number of sub-budgets

% Number of basis functions. This is used to calculate pi_hat.  Code
% currently uses polynomial basis functions (as opposed to splines).
poly_degree = polynomial_degree;  

% Number of bootstrap repetitions. 
% Time taken should be roughly linear in bootstrap_reps.
% New seed is set for every 1000 increase bootstrap_reps - so
% set bootstrap_reps to be the nearest thousandth. 
bootstrap_reps = bs_reps;      
bootstrap_reps = max(1000,1000*round(bootstrap_reps/1000));

% Set seed for reproducibility.  Need a new seed for each
% (budget_l,classes, poly_degree, bootstrap_reps, flag_estimator) case that
% one tries to ensure independence between results
seed = 10000*(bootstrap_reps/1000) + poly_degree*1000  + budget_l*100 + classes*10 + flag_estimator;

% Pre-specify number of cores to be used.  Truncate if this number exceeds 
% number of cores on the computer/server.
if feature('numCores') < cores
    num_cores = feature('numCores');
else
    num_cores = cores;
end

% Check to see if Parellel Computing installed
v_=ver;
[installedToolboxes{1:length(v_)}] = deal(v_.Name);
Parellel = all(ismember('Parallel Computing Toolbox',installedToolboxes));

% If installed, begin toolbox
if Parellel == 1
    delete(gcp);
    parpool(num_cores);
end

%% Create budgets 
% RUM_21_budgets cleans the dataset and outputs normalized budgets, as well
% as expenditure shares and income for each surveyed household. 
[budgets_all,shares_all,income_all,instrument_all] = RUM_21_budgets;


%% Re-map budgets_all to sub-budgets{budget_n} 
% budgets{i} is the ith sub-budget.
% I tried remapping shares_all and income_all, but ran into difficulty.
% This step does not really matter, it just makes the code a bit easier to
% follow.
budgets{budget_n} = [];
for ii = 1:budget_n
    budgets{ii} = budgets_all(ii:ii+(budget_l-1),:);
end

%% Compute X and A
% We compute X and A for each budget{ii} 
% Only do this if flag_genAX = 1.  Otherwise load X and A

    X{budget_n} = [];
    Tau_Set{budget_n} = [];
    fprintf('Starting to generate X and A using %d-length budgets \n',budget_l)
    for ii = 1:budget_n
        X{ii} = RUMCG_31_genX(budgets{ii});
        Tau_Set{ii} = RUMCG_PatternSet(tau_set_size,X{ii});
        fprintf('Completed budget %d out of %d. \n',ii,budget_n);
    end

Time_taken_AX = toc;


%% Pihat, Jstat, Critical Value, Pr(Jstat = 0)
% In this step we obtain a test statistic (Jstat) that is equal to zero
% iff Av = pi 
% iff we can find a distribution of pseudo-agents (v) that 
% satisfy such that the distribution matches the observed distribution of
% patch choices in the data (pi)
% iff household preferences are consistent with GARP+RUM
% iff Jstat = 0

% indJ =1 if we want to compute Jval and boostrap.
if indJ
pi_hat{budget_n} = [];
eta_hat{budget_n} = [];
Jstat_bs{budget_n} = [];
Jstat   = zeros(budget_n,1);
CV      = zeros(budget_n,2);
pval    = zeros(budget_n,1);
tau     = zeros(budget_n,1);

    % Compute Jstat and bootstrap distribution using Basis Smoothing
    for ii = 1:budget_n
        [pi_hat{ii},Jstat(ii,1),Jstat_bs{ii},CV(ii,:), pval(ii,1), tau(ii,1),eta_hat{ii}] ...
            =  RUMCG_71_EndogenousStatistic(budgets{ii},cell(shares_all(ii:ii+(budget_l-1))),cell(income_all(ii:ii+(budget_l-1))),cell(instrument_all(ii:ii+(budget_l-1))),X{ii},Tau_Set{ii});
        % Input of above function is budgets in period ii, the structure of
        % shares which is [num_household,classes]{budget_l} and income, which
        % is [num_household,1]{budget_l}.  Also input A and X for budget-period
        % ii.
        fprintf('Completed bootstrap %d out of %d. \n',ii,budget_n);
    end

end    

%% Save
Time_taken = toc;
c = clock;
datetime = strcat(num2str(c(1)),num2str(c(2)),num2str(c(3)),'_',num2str(c(4)),num2str(c(5)));
if flag_estimator == 1
     name = strcat('Output/RUM_',num2str(budget_l),'LengthBudget_',num2str(classes),'Classes_',num2str(poly_degree),'PolyDegree_',num2str(bootstrap_reps),'BSreps_Series_',num2str(tau_ind),'taunot0_',datetime);
elseif flag_estimator == 0
     name = strcat('Output/RUM_',num2str(budget_l),'LengthBudget_',num2str(classes),'Classes_',num2str(bootstrap_reps),'BSreps_Kernel_',num2str(tau_ind),'taunot0_',datetime);
else
     name = strcat('Output/RUM_',num2str(budget_l),'LengthBudget_',num2str(classes),'Classes_',num2str(poly_degree),'PolyDegree_',num2str(bootstrap_reps),'BSreps_Series_Endogenous',num2str(tau_ind),'taunot0_',datetime);
end
save(name); 



