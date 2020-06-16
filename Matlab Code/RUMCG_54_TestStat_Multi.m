function [JStat] = RUMX_54_TestStat_Multi(A, pi_hat,N,poly_degree,X,repetitions,TestStat,origin,Omega)
%% Setup
[I,H]=size(A);
J = size(X, 2);
first_patch = ones(size(X,2), 1);
for ii = 1:size(X,2)
    Nr_Patch(ii,1) = size(find(X(:,ii) == 0),1);
    first_patch(ii+1, 1) = first_patch(ii, 1) + Nr_Patch(ii, 1);
end
if nargin<9
    Omega=eye(I);
end
if nargin<8
    origin=zeros(I, 1);
end

% Check to see if Omega is PD and symmetric
[~,p] = chol(Omega);
if p > 0
    error('Omega is not positive definite')
end
if Omega ~= Omega.'
    error('Omega is not symmetric')
end

for ii = 1:repetitions
    pi_hat(:,ii) = pi_hat(:,ii) - origin;
end
% N = sum N_j/poly_degree
% We rescale N by poly_degree as per series estimation section in KS
% That is, we replace the vector N with N/max_jK(j), where K(j) is the
% polynomial degree in year j.  This needs to be adjusted if one uses a
% different polynomial degree for each year, in which case we use the
% largest polynomial degree.
N = sum(N)/poly_degree;

%% Build the master and pricing problems.
% Here we build the master problem.
Q = [Omega zeros(I, H)];
Q = [Q ; zeros(H, I + H)];
Amaster = [Omega A];
lb = [ones(I, 1)*-1; zeros(H, 1)];
f = zeros(H + I, 1);

% This matrix is used to check whether a given variable is actually used.
% for each master problem in which a variable is not used, it loses one
% credit. If it has zero credits remaining, it can be overwritten. 
% Credits = ones(H, 1)*50;

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


% This matrix is used to check whether a given variable is actually used.
% for each master problem in which a variable is not used, it loses one
% credit. If it has zero credits remaining, it can be overwritten. 
Credits = ones(H, 1)*50;

%% Main Algoritm
start_types = H;
added_types = 0;
JStat = zeros(repetitions, 1);
time_QP = 0;
time_IP = 0;
IP_Count = 0;
time_Heur = 0;
time_expand = 0;
for ii = 1:repetitions
fmaster = 999;
add_iterations = 0;
improve = 1;
while fmaster > 0 && improve == 1
    LB_Solution = 0;
    tic
    [xmaster, fmaster] = cplexqp(Q, f, [], [], Amaster, pi_hat(:,ii), lb, []);
    time_QP = toc + time_QP;
    distance = xmaster(1:I);
    nuhat = pi_hat(:,ii) - distance;
    JStat(ii)  =N*fmaster*2;
    %fprintf('Jstat(%d) = %d \n', ii, JStat(ii))
    if JStat(ii) < TestStat
        %fprintf('UB Break \n')
        break;
    end
    target = sum(distance.*nuhat);
    improve = 0;
    if fmaster > 0
        fprice = -99;
        add_iterations = add_iterations+1;
        if mod(add_iterations, 10) ~= 0
        tic
        [xprice, fprice] = RUMCG_55_PriceHeurBI(X, distance, target); 
        time_Heur = time_Heur + toc;
        end
        if fprice > target + 0.00001 && mod(add_iterations, 10) ~= 0 
            improve = 1;
        else
            priceobj = [-distance; zeros(J^2, 1)];
            tic
            [xprice, fprice, exitflag] = cplexbilp(priceobj, ConstraintIneq, RHSineq, ConstraintEq, RHSeq);
            time_IP = time_IP + toc;
            IP_Count = IP_Count +1;
            if -fprice > target +0.00001
            improve = 1;
            end
        end
    end
   
    if improve == 1
            % Every 100 times we add additional types, we expand the matrices. Expanding a matrix
            % takes time, but using matrices that are much larger than needed
            % would incur more overhead. 
            if mod(added_types, 100) == 0
                Amaster = [Amaster zeros(I, 100)];
                lb = [lb; zeros(100, 1)];
                f = [f; zeros(100, 1)];
                Q = [Q zeros(size(Q, 1), 100)];
                Q = [Q ; zeros(100,size(Q,2))];
                Credits = [Credits; zeros(100, 1)*50];
            end
            added_types = added_types + 1;
            Amaster(:,I+start_types+added_types) = xprice(1:I); 
    end
end
nuhat = nuhat + origin;
pi_hat(:,ii) = pi_hat(:,ii) + origin;
JStat(ii)  =N*(nuhat-pi_hat(:,ii)).'*Omega*(nuhat-pi_hat(:,ii));
if mod(ii, 100) == 0
    fprintf('Computed JStat %d \n', ii)
    end 
end

end