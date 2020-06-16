function X = RUMX_31_genX(budgets)
%% This function creates the Matrix X with a branching scheme.
% This algorithms works recursively.
% We test whether a patch can exist with the given relations to other budgets, 
% up to the given depth. (i.e. if the depth is 3, we will only look a the last 3 budgets).
% If a solution exists, we will continue by setting the budget just beyond the considered depth to -1 or 1. 
% This determines whether we test for the patch being above or below that budget.
% We then increase the depth and run Test_Patch again. 
global budget_l
X = zeros(0,0);
for ii = 1:budget_l
    Xrow = zeros(budget_l, 1);
   
    if ii == 1
        Xrow(2) = 1;
        X = RUMX_31_recurse(ii, budgets, X, Xrow, 2);
        Xrow(2) = -1;
        X = RUMX_31_recurse(ii, budgets, X, Xrow, 2);
    else
        Xrow(1) = 1;
        X = RUMX_31_recurse(ii, budgets, X, Xrow, 1);
        Xrow(1) = -1;
        X = RUMX_31_recurse(ii, budgets, X, Xrow, 1);
    end
end
end

function X = RUMX_31_recurse(current_budget, budgets, X, Xrow, depth);

global budget_l

    obj = budgets(current_budget, :)';
    Aineq = zeros(depth,size(budgets, 2));
    bineq = zeros(depth,1);
    lb = zeros(size(budgets, 2), 1);
    for ii = 1:depth
        if Xrow(ii) == 1
            Aineq(ii,:) = -budgets(ii,:);
            bineq(ii) = -1; 
        elseif Xrow(ii) == -1
            Aineq(ii,:) = budgets(ii,:);
            bineq(ii) = 1;
        end
    end
    [sol,~,exitflag] = cplexlp(obj, Aineq, bineq, obj', 1, lb, []);
    % The exitflag is > 0 if there is a feasible solution. 
    if depth == budget_l || (depth == budget_l -1 && current_budget == budget_l)
        if exitflag > 0
            X = [X ; Xrow'];
        end
    else
        if exitflag > 0
            if depth == current_budget - 1
                Xrow(depth + 2) = 1;
                X = RUMX_31_recurse(current_budget, budgets, X, Xrow, depth+2);
                Xrow(depth + 2) = -1;
                X = RUMX_31_recurse(current_budget, budgets, X, Xrow, depth+2);
            else
                Xrow(depth + 1) = 1;
                X = RUMX_31_recurse(current_budget, budgets, X, Xrow, depth+1);
                Xrow(depth + 1) = -1;
                X = RUMX_31_recurse(current_budget, budgets, X, Xrow, depth+1);
            end
        end
    end
end