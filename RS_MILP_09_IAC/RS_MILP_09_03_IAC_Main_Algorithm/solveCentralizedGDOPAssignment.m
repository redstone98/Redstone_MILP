function result = solveCentralizedGDOPAssignment(time_index, ...
                                                topK_GDOP_history, ...
                                                topK_NRHO_ID_history, ...
                                                topK_L4_ID_history, ...
                                                topK_L5_ID_history, ...
                                                capacity_NRHO, ...
                                                capacity_L4, ...
                                                capacity_L5)

%% ========================================================================
% Centralized Top-K GDOP Assignment
% ========================================================================

N_FRT = size(topK_GDOP_history,2);
K     = size(topK_GDOP_history,3);

N_NRHO = length(capacity_NRHO);
N_L4   = length(capacity_L4);
N_L5   = length(capacity_L5);


%% ------------------------------------------------------------------------
% Extract timestep
% -------------------------------------------------------------------------

G = reshape( ...
    topK_GDOP_history(time_index,:,:), ...
    N_FRT,K);

NRHO_ID = reshape( ...
    topK_NRHO_ID_history(time_index,:,:), ...
    N_FRT,K);

L4_ID = reshape( ...
    topK_L4_ID_history(time_index,:,:), ...
    N_FRT,K);

L5_ID = reshape( ...
    topK_L5_ID_history(time_index,:,:), ...
    N_FRT,K);


%% ------------------------------------------------------------------------
% Variables
%
% x(i,k)
%
% Linear index:
% (i-1)*K + k
% -------------------------------------------------------------------------

nVar = ...
    N_FRT * K;

f = zeros(nVar,1);

lb = zeros(nVar,1);
ub = ones(nVar,1);


for i = 1:N_FRT

    for k = 1:K

        idx = ...
            (i-1)*K + k;

        if isfinite(G(i,k))

            f(idx) = ...
                G(i,k);

        else

            % Invalid candidate
            f(idx) = 0;
            ub(idx) = 0;

        end

    end
end


%% ========================================================================
% Each FRT selects exactly one candidate
% ========================================================================

Aeq = zeros(N_FRT,nVar);

beq = ones(N_FRT,1);


for i = 1:N_FRT

    idx = ...
        (i-1)*K + (1:K);

    Aeq(i,idx) = 1;

end


%% ========================================================================
% Servicing satellite capacity constraints
% ========================================================================

nConstraint = ...
    N_NRHO + ...
    N_L4 + ...
    N_L5;

A = zeros(nConstraint,nVar);

b = [
    capacity_NRHO(:)
    capacity_L4(:)
    capacity_L5(:)
    ];


row = 0;


%% NRHO capacity

for s = 1:N_NRHO

    row = row + 1;

    for i = 1:N_FRT
        for k = 1:K

            if NRHO_ID(i,k) == s

                idx = ...
                    (i-1)*K + k;

                A(row,idx) = 1;

            end
        end
    end
end


%% L4 capacity

for s = 1:N_L4

    row = row + 1;

    for i = 1:N_FRT
        for k = 1:K

            if L4_ID(i,k) == s

                idx = ...
                    (i-1)*K + k;

                A(row,idx) = 1;

            end
        end
    end
end


%% L5 capacity

for s = 1:N_L5

    row = row + 1;

    for i = 1:N_FRT
        for k = 1:K

            if L5_ID(i,k) == s

                idx = ...
                    (i-1)*K + k;

                A(row,idx) = 1;

            end
        end
    end
end


%% ========================================================================
% Solve MILP
% ========================================================================

intcon = ...
    1:nVar;

options = optimoptions( ...
    'intlinprog', ...
    'Display','off');


[x,fval,exitflag,output] = ...
    intlinprog( ...
        f, ...
        intcon, ...
        A,b, ...
        Aeq,beq, ...
        lb,ub, ...
        options);


%% ========================================================================
% Decode solution
% ========================================================================

result.success = ...
    exitflag > 0;

result.objective = ...
    fval;

result.exitflag = ...
    exitflag;

result.output = ...
    output;


result.rank = NaN(N_FRT,1);
result.GDOP = NaN(N_FRT,1);

result.NRHO_ID = NaN(N_FRT,1);
result.L4_ID   = NaN(N_FRT,1);
result.L5_ID   = NaN(N_FRT,1);


if result.success

    for i = 1:N_FRT

        idx = ...
            (i-1)*K + (1:K);

        selected = ...
            find(x(idx) > 0.5,1);

        result.rank(i) = ...
            selected;

        result.GDOP(i) = ...
            G(i,selected);

        result.NRHO_ID(i) = ...
            NRHO_ID(i,selected);

        result.L4_ID(i) = ...
            L4_ID(i,selected);

        result.L5_ID(i) = ...
            L5_ID(i,selected);

    end
end

end