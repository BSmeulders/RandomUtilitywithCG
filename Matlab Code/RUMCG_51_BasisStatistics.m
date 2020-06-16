function [pi_hat, Jstat, Jstat_bs, CV, prob,tau,nuhat] =  RUMX_51_BasisStatistics(budgets,shares,income,X,A) 
%% Code Description: Test statistic for H_0: Households satisfy GARP + RUM
% This module computes pihat (probability of choosing patch i from budget j),
% the Jstat, and the boostrapped Jstat distrubiton
% using method described in Bootstrap Algorithm with Tightening
% First, compute pihat using method described in series estimation.
% Second, get untightened Jstat and tightened nu_hat using tau.
% Third, re-draw the data and calculate a new bootstrapped pihat.  
% Recenter the tightened nuhat using (pihat_bs - pihat).  
% Last, get critical values and probabilities for the untightened Jstat
% Input:
%   - budgets: median-normalized budget constraints
%   - shares:  Expenditure share for each person n
%   - income:  Total expenditure on non-durables for each person n
%   - X:       Matrix defining patches
%   - A:       Matrix of GARP consistent pseudo-agents 
% Output:
%   - pi_hat:   Probability of purchasing from patch i on budget j
%   - Jstat:    The tau-tightened J statistic
%   - Jstat_bs: Bootstrapped distribution of Jstat
%   - CV:       CV(1) is the 95 percentile critical value and CV(2) is the
%               99 percentile
%  - prob:      Probability that the tau-tightened J statistic is zero
%  - tau:       Tuning parameter 

%% Global parameters
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores


%% Log income
for jj = 1:budget_l
   income{jj} = log(income{jj}); 
end

%% Median income
% Calculate median income.  This is held fixed through the bootstrap, so it
% is necessary to do it ex ante.
% I let med_income be a structure rather than a matrix, because (needs to
% be adjusted in the code) we can potentially allow for a  non-constant
% polynomial degree over the different years.
med_income{budget_l} = [];    % Polynomial basis applied to medium income
for jj = 1:budget_l
    med_income{jj} = 1;
    for pp = 1:poly_degree
        med_income{jj} = [med_income{jj}, median(income{jj}).^pp];
    end
end

%% Pihat
% We calculate pi_hat, which is the probability of a household with
% median income purchasing patch i from budget j.
[pi_hat,N] =RUM_53_BasisPihat(budgets,shares,income,med_income,X,budget_l,poly_degree);

%% Tau: The Tuning parameter
% To get tau, we first need to get an estimate of the variance v_N.
% We get this via bootstrapping the variance of pi_{EC}.  
% Code for bootstrap is a little complicated, because:
%    - Need to preset cores and check if parellel computing is installed
%       to use ECCO
%    - Need to avoid slicing
%    - Need to set seeds correctly (hence need to use subseeds)
%    - MatLab strugles with structures in parfor loops, so we need to 
%       randomize in a subfunction
% Preassign varaibles for bootstraps
J = ceil(bootstrap_reps/num_cores);
I = num_cores;
B = bootstrap_reps;
size_pi = size(pi_hat,1);
pi_bs_temp{I} = [];
pi_bs = zeros(size(pi_hat,1),J*I);
Jstat_bs_temp{I} = [];
Jstat_bs  = zeros(J*I,1);

parfor ii = 1:I
    pi_bs_temp2 = zeros(size_pi,J);
    for jj = 1:J
        % Effective bootstrap iteration
        bb = (ii-1)*J + jj

        % If we have reached total bootstrap_reps, continue
        if bb > B
            continue;
        end
        
        % Set seed: Since we are doing this in parallel we need to set seed on 
        % each loop.  To ensure independence we set subseed to be bb rather 
        % than the seed.
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = bb;

        % Resample shares and income
        % We need to run this sub-function because MatLab cannot handle
        % structures in parallel.
        [shares_bs,income_bs] = RUM_52_BasisRandomize(shares,income,N,budget_l);

        % Get new pi_hat based on resampled shares and income, but do not 
        % update median income.
        pi_bs_temp2(:,jj) =RUM_53_BasisPihat(budgets,shares_bs,income_bs,med_income,X,budget_l,poly_degree);   
    end
    pi_bs_temp{ii} = pi_bs_temp2;
end
for ii=1:I
    pi_bs(:, (ii-1)*J+1:ii*J) = pi_bs_temp{ii};
end
if B ~= I*J
    pi_bs(:,B+1:I*J) = [];
end

% For each pi_{i,j} calculate variance.  We do not need to worry about
% covariance.
var_bs = var(pi_bs.',1).';

% For each year jj, estimate of variance v_N and hence 
% n_K = min_j N_j I_j/trace(v_Nj)
n_K = zeros(budget_l,1);
for jj = 1:budget_l
   v_N =  N(jj)*var_bs(find(X(:,jj) == 0));
   I_j = size(find(X(:,jj) == 0),1);
   n_K(jj,1) = N(jj)*I_j/sum(v_N);
end

% n_K is the minimum n_K over the budget years.
n_K = min(n_K);

% Tau: The Tuning parameter
tau = (log(n_K)/n_K)^0.5; % Tau is fixed for bootstrap

% Given Tau, we can compute the origin of the cone that is used for
% computing the test statistic distribution.
origin = sum(A, 2)/size(A,2)*tau;
%% Jstat and bootstrap distribution
% Jstat is not tightened
[nuhat,Jstat,A] =  RUMX_54_BasisTeststat(A,pi_hat,N,poly_degree,X);

% Get tau-tightened nuhat
[nuhat_tight,~,A] =  RUMX_54_BasisTeststat(A,pi_hat,N,poly_degree,X,origin);

% Resample shares income to get non-parametric pihat. 
Jstat_bs = zeros(B,1);
for ii = 1:B
    % Resample shares and income
    [shares_bs,income_bs] = RUM_52_BasisRandomize(shares,income,N,budget_l);
    % Get new pi_hat based on resampled shares and income, but do not 
    % update median income.
    pi_bs2 =RUM_53_BasisPihat(budgets,shares_bs,income_bs,med_income,X,budget_l,poly_degree); 
    % Recenter the tau-tightned nuhat
    pi_bs2 = pi_bs2 - pi_hat + nuhat_tight;

    % Get the Jstatistic associated with the recenter the tau-tightned nuhat
    % Do not update tau.
    [~,Jstat_bs(ii,1),A] =  RUMX_54_BasisTeststat(A,pi_bs2,N,poly_degree,X,origin);
    fprintf('Iteration %d \n', ii)
end


%% Critical value and Pr(Jstat == 0)
Jstat_bs = sortrows(Jstat_bs,1);
cv_95 = Jstat_bs( ceil(bootstrap_reps*0.95));
cv_99 = Jstat_bs( ceil(bootstrap_reps*0.99));
CV = [cv_95  cv_99];
if isempty(min(Jstat_bs( Jstat_bs>Jstat)))
    prob = 0;
else
    prob  = 1- find(min(Jstat_bs( Jstat_bs>Jstat)) == Jstat_bs)/bootstrap_reps;
end



end