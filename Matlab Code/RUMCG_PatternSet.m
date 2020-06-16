function A = RUMX_PatternSet( nr_start_types, X )
%RUMX_PATTERNSET
% We generate a subset of the rational choice types. We do so in a
% semi-random way.

% Global Variables
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores

% We first compute how many patches there are on each budget.
Nr_Patch = zeros(budget_l,1);
for ii = 1:budget_l
    Nr_Patch(ii,1) = size(find(X(:,ii) == 0),1);
end
A = zeros(sum(Nr_Patch), nr_start_types);
for ii = 1:nr_start_types
    A(:,ii) = RUMX_PatternGeneration(X, Nr_Patch);  
end

end

function Aout = RUMX_PatternGeneration(X, Nr_Patch)
% Patterns are generated as follows. First, a random patch is chosen on
% each budget. Next, we test whether the resulting type is rational. If
% not, a minimal change is made to the type. On one budget per Strongly
% connected component, one new patch
% is chosen, that minimizes the number of changes, while removing that
% budget from a Strongly Connected Component. This is repeated until a
% rational choice type is found.

global budget_l

% Build an initial type.
Aout = zeros(sum(Nr_Patch),1);
sum_patch = 0;
for ii = 1:size(Nr_Patch)
    patch = sum_patch + ceil(rand(1)*Nr_Patch(ii));
    Aout(patch) = 1;
    sum_patch = sum_patch + Nr_Patch(ii);
end

valid_type = 0;
iterations = 0;
while valid_type == 0 && iterations < 1000
    
    % Generate a Revealed Preference Graph
    RPGraph = CreateRPGraph(X, Nr_Patch, Aout);
    % Compute the strongly connected components.
    [SCC, SCCsize, nrSCC] = Tarjans(RPGraph);

    % If there are SCC involving more than one patch, choose one buget
    % involved and change it (minimally) to remove RP relations.
    if nrSCC == budget_l
        valid_type = 1;
    end
    if nrSCC < budget_l
        for ii = 1:nrSCC
            if SCCsize(ii) > 1 % If there is an SSC with multiple budgets, we change the patch for one of these.
                change_budget = SCC(ii,ceil(rand(1)*SCCsize(ii)));
                % Compute what the first patch of the changed budget is.
                first_budget_patch = 1;
                for jj = 1:change_budget-1
                    first_budget_patch = first_budget_patch + Nr_Patch(jj);
                end
                % Save the number of the currently used patch for this
                % budget.
                for jj = first_budget_patch:first_budget_patch+Nr_Patch(change_budget)-1
                    if Aout(jj) == 1
                        current_patch = jj;
                    end;
                end
                patch_scores = zeros(Nr_Patch(change_budget), 1);
                for jj = 1:Nr_Patch(change_budget)
                    % First check whether the new patch is the same within
                    % the SCC.
                    no_pref_remove = 1;
                    for kk = 1:SCCsize(ii)
                        if X(current_patch, SCC(ii,kk)) < X(jj+first_budget_patch-1,SCC(ii,kk))
                            no_pref_remove = 0;
                        end
                    end
                    if no_pref_remove == 1
                        patch_scores(jj) = 999;
                    % If it isn't give it a score which is smaller the
                    % closer the new patch is to the current one. Removed
                    % RP relations cost less than added relations.
                    else
                        for kk = 1:budget_l
                            if X(current_patch,kk) == -1 && X(jj+first_budget_patch-1,kk) == 1
                                patch_scores(jj) = patch_scores(jj) + 1;
                            end
                            if X(current_patch,kk) == 1 && X(jj+first_budget_patch-1,kk) == -1
                                patch_scores(jj) = patch_scores(jj) + 5;
                            end
                        end
                    end
                end
                % Find the patch with the lowest score
                [min_score, min_patch] = min(patch_scores);
                Aout(current_patch) = 0;
                Aout(first_budget_patch + min_patch -1) = 1;
            end
        end
    end   
    iterations = iterations + 1;
end
% Find one patch involved in a violation and minimally change it (move to
% a patch that is as close as possible, but removes some RP relations).

end

function RPGraph = CreateRPGraph(X, Nr_Patch, Aout)

global budget_l

RPGraph = zeros(budget_l, budget_l);
budget = 0;
for ii = 1:size(Aout)
    if Aout(ii) == 1
        budget = budget + 1;
        for jj = 1:budget_l
          if X(ii, jj) == -1
            RPGraph(jj, budget) = 1;
          end
        end 
    end
end

end

function [SCC, SCCsize, nrSCC] = Tarjans(RPGraph)
global budget_l
index = zeros(budget_l, 1);
totindex = 1;
lowlink = zeros(budget_l, 1);
stack = zeros(budget_l, 1);
onstack = zeros(budget_l, 1);
stack_length = 0;
SCCsize = zeros(budget_l, 1);
SCC = zeros(budget_l, budget_l);
nrSCC = 0;

for ii = 1:budget_l
    if index(ii) == 0
    [index,totindex, lowlink, stack, onstack, stack_length, SCC, SCCsize, nrSCC] = strongconnect(index, totindex, lowlink, stack, onstack, stack_length, SCC, SCCsize, nrSCC, RPGraph, ii);
    end
end
end

function [index,totindex, lowlink, stack, onstack, stack_length, SCC, SCCsize, nrSCC] = strongconnect(index, totindex, lowlink, stack, onstack, stack_length, SCC, SCCsize, nrSCC, RPGraph, ii)
global budget_l

index(ii) = totindex; lowlink(ii) = totindex;
totindex = totindex + 1;
stack_length = stack_length + 1;
stack(stack_length) = ii;
onstack(ii) = 1;

for jj = 1:budget_l
    if RPGraph(ii,jj) == 1 && index(jj) == 0
        [index,totindex, lowlink, stack, onstack, stack_length, SCC, SCCsize, nrSCC] = strongconnect(index, totindex, lowlink, stack,onstack, stack_length, SCC, SCCsize,nrSCC, RPGraph, jj);
        lowlink(ii) = min(lowlink(ii), lowlink(jj));
    elseif RPGraph(ii,jj) == 1 && onstack(jj) == 1
        lowlink(ii) = min(lowlink(ii), index(jj));
    end
end

if lowlink(ii) == index(ii)
    nrSCC = nrSCC + 1;
    while ii ~= stack(stack_length)
    SCCsize(nrSCC) = SCCsize(nrSCC) + 1;
    SCC(nrSCC,SCCsize(nrSCC)) = stack(stack_length);
    onstack(stack(stack_length)) = 0;
    stack(stack_length) = 0;
    stack_length = stack_length - 1;
    end
    SCCsize(nrSCC) = SCCsize(nrSCC) + 1;
    SCC(nrSCC,SCCsize(nrSCC)) = stack(stack_length);
    onstack(stack(stack_length)) = 0;
    stack(stack_length) = 0;
    stack_length = stack_length - 1;
    
end

end
