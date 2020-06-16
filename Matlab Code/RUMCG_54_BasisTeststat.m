function [nuhat,fval, A] =  RUMX_54_BasisTeststat(A,pi_hat,N,poly_degree,X,origin,Omega) 
%% Code Description
% nuhat solves min (A*nu - pihat)'Omega(A*nu - pihat) subject to nu >= 0. 
% Input:
%   - A is the agent-type matrix that determines which patch a particular
%       type of agent picks from budget j.
%   - pihat is computed using data on patch choices
%   - N is a Jx1 vector with the number of people sampled from budget
%       j=1,...,J
%   - tau is the tuning parameter.  If not specified, it is set to 0.
%   -Omega is a consistent estimator for the asymptotic variance.  Default
%   is to set it equal to I.

%% Global parameters
% This function is passed through a parfor loop.  Therefore, we cannot
% include global variables.  Instead, I put the variables that we need as
% an argument.  

%% Setup
[I,H]=size(A);
J = size(X, 2);
for ii = 1:J
    Nr_Patch(ii,1) = size(find(X(:,ii) == 0),1);
end
if nargin<7
    Omega=eye(I);
end
if nargin<6
    origin=zeros(I, 1);
end

pi_hat = pi_hat - origin;
% Check to see if Omega is PD and symmetric
[~,p] = chol(Omega);
if p > 0
    error('Omega is not positive definite')
end
if Omega ~= Omega.'
    error('Omega is not symmetric')
end
% N = sum N_j/poly_degree
% We rescale N by poly_degree as per series estimation section in KS
% That is, we replace the vactor N with N/max_jK(j), where K(j) is the
% polynomial degree in year j.  This needs to be adjusted if one uses a
% different polynomial degree for each year, in which case we use the
% largest polynomial degree.
N = sum(N)/poly_degree;

%% Find nuhat
% Nuhat is computed by a column generation algorithm. In a first iteration, only the
% initially computed types are taken into account. Next, a pricing problem
% is solved, to check whether there exist additional types that can, if
% added to the problem, improve the solution. This is done repeatedly,
% until either it is proven that pi_hat lies inside the cone, or no more
% types are found to improve the solution.

% Here we build the master problem.
Q = [Omega zeros(I, H)];
Q = [Q ; zeros(H, I + H)];
Amaster = [Omega A];
lb = [ones(I, 1)*-1; zeros(H, 1)];
f = zeros(H + I, 1);

% Next we build the pricing problem. We have four kinds of constraints we
% begin with the equality constraints.
% 1) Assignment Constraint: When building a new type, we need at least one
% patch chosen on each budget.
ConstraintEq = zeros(J, I + J^2);
count = 1;
for ii = 1:J
    for jj = 1:Nr_Patch(ii);
        ConstraintEq(ii,count) = 1;
        count = count + 1;
    end
end
% 2) Anti-symmetry constraint. We enforce that either the patch on budget i
% is preferred over the patch on j or vice versa.
ConstraintAntiSymmetry = zeros(J*(J-1)/2, I + J^2);
count = 1;
for ii = 1:(J-1)
    for jj = ii+1:J
        ConstraintAntiSymmetry(count, I + (ii-1)*J + jj) = 1;
        ConstraintAntiSymmetry(count, I + (jj-1)*J + ii) = 1;
        count = count + 1;
    end
end
% 3) Revealed Preference Constraint: The choice of certain patches implies
% revealed preference constraints. 
ConstraintRP = zeros(J^2, I + J^2);
count = 1;
for ii = 1:J
    for jj = 1:Nr_Patch(ii)
        for kk = 1:J
            if X(count, kk) == -1
                ConstraintRP((kk-1)*J + ii, count) = 1;
            end
        end
        count = count + 1;
    end
end
for ii = 1:J^2
    ConstraintRP(ii, I + ii) = -1;
end
% 4) Transitivity constraints. If i preferred to j and j to k, then i must also be preferred to k.
ConstraintTrans = zeros(J*(J-1)*(J-2), I + J^2);
count = 1;
for ii = 1:J
    for jj = 1:J
        if ii ~= jj
            for kk = 1:J
                if ii ~= kk && jj ~= kk
                    ConstraintTrans(count, I + (ii-1)*J + jj) = 1;
                    ConstraintTrans(count, I + (jj-1)*J + kk) = 1;
                    ConstraintTrans(count, I + (ii-1)*J + kk) = -1;
                    count = count + 1;
                end
            end
        end
    end
end
% Merge the inequality constraints.
ConstraintIneq = [ConstraintAntiSymmetry ; ConstraintRP; ConstraintTrans];

% Set the right hand sides of the Constraints
RHSeq = ones(J, 1);
RHSineq = [ones(J*(J-1)/2, 1); zeros(J^2, 1); ones(J*(J-1)*(J-2), 1)];


%% Column Generation
fmaster = 999;
improve = 1;
iterations = 0;
while fmaster > 0 && improve == 1
    [xmaster, fmaster] = cplexqp(Q, f, [], [], Amaster, pi_hat, lb, []);  
    distance = xmaster(1:I);
    nuhat = pi_hat - distance;
    target = sum(distance.*nuhat);
    improve = 0;
    if fmaster > 0
        [xprice, fprice] = RUMCG_55_PriceHeurBI(X, distance, target);
        if fprice > target + 0.00001 
            improve = 1;
        else
            priceobj = [-distance; zeros(J^2, 1)];
            [xprice, fprice, exitflag] = cplexbilp(priceobj, ConstraintIneq, RHSineq, ConstraintEq, RHSeq);
            if -fprice > target +0.00001
            improve = 1;
            end
        end
    end
    if improve == 1
        % Every 100 iterations, we expand the matrices. Expanding a matrix
        % takes time, but using matrices that are much larger than needed
        % would incur more overhead. 
        if mod(iterations, 100) == 0
            fprintf('Iteration %d: %d \n',iterations, fmaster)
            Amaster = [Amaster zeros(I, 100)];
            lb = [lb; zeros(100, 1)];
            f = [f; zeros(100, 1)];
            Q = [Q zeros(size(Q, 1), 100)];
            Q = [Q ; zeros(100,size(Q,2))];
        end
        iterations = iterations + 1;
        Amaster(:,I+size(A,2)+iterations) = xprice(1:I);      
    end
end

%% Output
% We save the columns added during the columns generation, for use later
% on.
A = Amaster(:,I+1:I+size(A,2)+iterations);
nuhat = nuhat + origin;
pi_hat = pi_hat + origin;
fval  =N*(nuhat-pi_hat).'*Omega*(nuhat-pi_hat);

end