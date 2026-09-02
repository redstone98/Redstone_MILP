function result = solveDistributedSequential_BCD( ...
    time_index, ...
    topK_GDOP_history, ...
    topK_NRHO_ID_history, ...
    topK_L4_ID_history, ...
    topK_L5_ID_history, ...
    capacity_NRHO, ...
    capacity_L4, ...
    capacity_L5, ...
    decision_order)

%% ========================================================================
% Sequential Distributed Local-MILP Assignment
%
% Cyclic best-response / block-coordinate integer optimization
%
% At each FRT update:
%
%   1. Fix all other FRT decisions
%   2. Release current FRT resources
%   3. Solve local integer problem for current FRT
%   4. Update its decision
%   5. Hand over to next FRT
%
% Convergence:
%
%   Stop when one complete sweep produces no decision changes.
%
% Optional:
%
%   initial_rank : N_FRT x 1 initial feasible assignment
%   max_sweeps   : maximum number of complete sweeps
%
% Old calling format remains valid:
%
% result = solveDistributedSequential( ...
%     time_index, ...
%     topK_GDOP_history, ...
%     topK_NRHO_ID_history, ...
%     topK_L4_ID_history, ...
%     topK_L5_ID_history, ...
%     capacity_NRHO, ...
%     capacity_L4, ...
%     capacity_L5, ...
%     decision_order)
%
% ========================================================================

%% ------------------------------------------------------------------------
% Problem size
% -------------------------------------------------------------------------

N_FRT = size(topK_GDOP_history,2);
K     = size(topK_GDOP_history,3);

N_NRHO = length(capacity_NRHO);
N_L4   = length(capacity_L4);
N_L5   = length(capacity_L5);

capacity_NRHO = capacity_NRHO(:);
capacity_L4   = capacity_L4(:);
capacity_L5   = capacity_L5(:);

decision_order = decision_order(:);


%% ------------------------------------------------------------------------
% Optional inputs
% -------------------------------------------------------------------------

if nargin < 10
    initial_rank = [];
end

if nargin < 11 || isempty(max_sweeps)
    max_sweeps = 100;
end


%% ------------------------------------------------------------------------
% Check decision order
% -------------------------------------------------------------------------

if length(decision_order) ~= N_FRT || ...
        ~isequal(sort(decision_order),(1:N_FRT)')

    error(['decision_order must contain each FRT index ' ...
           'exactly once.']);

end


%% ------------------------------------------------------------------------
% Extract current timestep
% -------------------------------------------------------------------------

G = reshape( ...
    topK_GDOP_history(time_index,:,:), ...
    N_FRT,K);

NRHO = reshape( ...
    topK_NRHO_ID_history(time_index,:,:), ...
    N_FRT,K);

L4 = reshape( ...
    topK_L4_ID_history(time_index,:,:), ...
    N_FRT,K);

L5 = reshape( ...
    topK_L5_ID_history(time_index,:,:), ...
    N_FRT,K);


%% ------------------------------------------------------------------------
% intlinprog options
% -------------------------------------------------------------------------

options = optimoptions( ...
    'intlinprog', ...
    'Display','off');


%% ========================================================================
% Initialization
% ========================================================================

load_NRHO = zeros(N_NRHO,1);
load_L4   = zeros(N_L4,1);
load_L5   = zeros(N_L5,1);

selected_rank = NaN(N_FRT,1);


%% ------------------------------------------------------------------------
% Case 1:
% User supplies an initial feasible assignment
% -------------------------------------------------------------------------

if ~isempty(initial_rank)

    initial_rank = initial_rank(:);

    if length(initial_rank) ~= N_FRT

        error('initial_rank must have length N_FRT.');

    end


    for i = 1:N_FRT

        k = initial_rank(i);

        if ~isfinite(k) || ...
                k < 1 || ...
                k > K || ...
                k ~= round(k)

            error('Invalid initial_rank for FRT %d.',i);

        end


        if ~isfinite(G(i,k))

            error( ...
                'Initial candidate for FRT %d is invalid.', ...
                i);

        end


        n  = NRHO(i,k);
        l4 = L4(i,k);
        l5 = L5(i,k);


        load_NRHO(n) = load_NRHO(n) + 1;
        load_L4(l4)  = load_L4(l4) + 1;
        load_L5(l5)  = load_L5(l5) + 1;

        selected_rank(i) = k;

    end


    if any(load_NRHO > capacity_NRHO) || ...
       any(load_L4   > capacity_L4)   || ...
       any(load_L5   > capacity_L5)

        error('initial_rank is not globally feasible.');

    end


%% ------------------------------------------------------------------------
% Case 2:
% No initial assignment
%
% Generate a feasible initial point sequentially.
% Each FRT solves its own local MILP.
% -------------------------------------------------------------------------

else

    for order_index = 1:N_FRT

        i = decision_order(order_index);


        [local_success, k_new, ~, ~] = ...
            solveLocalFRTMILP( ...
                i, ...
                G, ...
                NRHO, ...
                L4, ...
                L5, ...
                load_NRHO, ...
                load_L4, ...
                load_L5, ...
                capacity_NRHO, ...
                capacity_L4, ...
                capacity_L5, ...
                options);


        if ~local_success

            result.success = false;
            result.converged = false;

            result.status = ...
                sprintf( ...
                'Initialization failed at FRT %d.', ...
                i);

            result.rank = selected_rank;

            result.GDOP = NaN(N_FRT,1);

            result.NRHO_ID = NaN(N_FRT,1);
            result.L4_ID   = NaN(N_FRT,1);
            result.L5_ID   = NaN(N_FRT,1);

            result.objective = NaN;

            result.load_NRHO = load_NRHO;
            result.load_L4   = load_L4;
            result.load_L5   = load_L5;

            return;

        end


        selected_rank(i) = k_new;


        n  = NRHO(i,k_new);
        l4 = L4(i,k_new);
        l5 = L5(i,k_new);


        load_NRHO(n) = ...
            load_NRHO(n) + 1;

        load_L4(l4) = ...
            load_L4(l4) + 1;

        load_L5(l5) = ...
            load_L5(l5) + 1;

    end

end


%% ========================================================================
% Initial objective
% ========================================================================

objective_history = NaN(max_sweeps+1,1);

objective_history(1) = ...
    computeAssignmentObjective( ...
        G, ...
        selected_rank);


rank_history = ...
    NaN(N_FRT,max_sweeps+1);

rank_history(:,1) = ...
    selected_rank;


%% ========================================================================
% Cyclic Sequential Best Response
% ========================================================================

converged = false;

local_solve_count = 0;

actual_sweeps = 0;


for sweep = 1:max_sweeps

    actual_sweeps = sweep;


    %% ---------------------------------------------------------------
    % Save decisions before complete sweep
    % ---------------------------------------------------------------

    previous_rank = ...
        selected_rank;


    %% ---------------------------------------------------------------
    % Sequential handover
    % ---------------------------------------------------------------

    for order_index = 1:N_FRT

        i = ...
            decision_order(order_index);


        %% ------------------------------------------------------------
        % Current assignment
        % ------------------------------------------------------------

        k_old = ...
            selected_rank(i);

        n_old  = NRHO(i,k_old);
        l4_old = L4(i,k_old);
        l5_old = L5(i,k_old);


        %% ------------------------------------------------------------
        % Release current FRT resources
        %
        % After this operation, load variables represent ONLY
        % decisions of all other FRTs.
        % ------------------------------------------------------------

        load_NRHO(n_old) = ...
            load_NRHO(n_old) - 1;

        load_L4(l4_old) = ...
            load_L4(l4_old) - 1;

        load_L5(l5_old) = ...
            load_L5(l5_old) - 1;


        %% ------------------------------------------------------------
        % Solve local integer problem
        %
        % All other FRT decisions are fixed.
        % ------------------------------------------------------------

        [ ...
            local_success, ...
            k_candidate, ...
            local_objective, ...
            exitflag] = ...
            solveLocalFRTMILP( ...
                i, ...
                G, ...
                NRHO, ...
                L4, ...
                L5, ...
                load_NRHO, ...
                load_L4, ...
                load_L5, ...
                capacity_NRHO, ...
                capacity_L4, ...
                capacity_L5, ...
                options);


        local_solve_count = ...
            local_solve_count + 1;


        %% ------------------------------------------------------------
        % Local problem should remain feasible because the old
        % assignment is feasible after releasing itself.
        % ------------------------------------------------------------

        if ~local_success

            % Restore previous assignment
            load_NRHO(n_old) = ...
                load_NRHO(n_old) + 1;

            load_L4(l4_old) = ...
                load_L4(l4_old) + 1;

            load_L5(l5_old) = ...
                load_L5(l5_old) + 1;


            result.success = false;
            result.converged = false;

            result.status = ...
                sprintf( ...
                ['Local MILP failed for FRT %d ' ...
                 'at sweep %d, exitflag = %d.'], ...
                i, ...
                sweep, ...
                exitflag);

            result.rank = selected_rank;

            return;

        end


        %% ------------------------------------------------------------
        % Tie handling
        %
        % If the old decision has the same GDOP as the local optimum,
        % retain the old decision.
        %
        % This prevents unnecessary switching/cycling between
        % equal-cost candidates.
        % ------------------------------------------------------------

        old_objective = ...
            G(i,k_old);

        comparison_tolerance = ...
            1e-10 * ...
            max([ ...
                1, ...
                abs(old_objective), ...
                abs(local_objective)]);


        if old_objective <= ...
                local_objective + comparison_tolerance

            k_new = ...
                k_old;

        else

            k_new = ...
                k_candidate;

        end


        %% ------------------------------------------------------------
        % Update current FRT decision
        % ------------------------------------------------------------

        selected_rank(i) = ...
            k_new;


        n_new  = NRHO(i,k_new);
        l4_new = L4(i,k_new);
        l5_new = L5(i,k_new);


        load_NRHO(n_new) = ...
            load_NRHO(n_new) + 1;

        load_L4(l4_new) = ...
            load_L4(l4_new) + 1;

        load_L5(l5_new) = ...
            load_L5(l5_new) + 1;


        % Next FRT now sees this updated assignment.
        %
        % This is a Gauss-Seidel / sequential handover update.

    end


    %% ---------------------------------------------------------------
    % End of one complete sweep
    % ---------------------------------------------------------------

    objective_history(sweep+1) = ...
        computeAssignmentObjective( ...
            G, ...
            selected_rank);

    rank_history(:,sweep+1) = ...
        selected_rank;


    %% ---------------------------------------------------------------
    % Convergence condition
    %
    % Stop ONLY after a complete sweep if every FRT decision
    % is identical to the previous sweep.
    % ---------------------------------------------------------------

    if isequal( ...
            selected_rank, ...
            previous_rank)

        converged = true;
        break;

    end

end


%% ========================================================================
% Trim history
% ========================================================================

objective_history = ...
    objective_history(1:actual_sweeps+1);

rank_history = ...
    rank_history(:,1:actual_sweeps+1);


%% ========================================================================
% Decode final solution
% ========================================================================

selected_GDOP = NaN(N_FRT,1);

selected_NRHO = NaN(N_FRT,1);
selected_L4   = NaN(N_FRT,1);
selected_L5   = NaN(N_FRT,1);


for i = 1:N_FRT

    k = ...
        selected_rank(i);

    selected_GDOP(i) = ...
        G(i,k);

    selected_NRHO(i) = ...
        NRHO(i,k);

    selected_L4(i) = ...
        L4(i,k);

    selected_L5(i) = ...
        L5(i,k);

end


%% ========================================================================
% Final feasibility check
% ========================================================================

capacity_feasible = ...
    all(load_NRHO <= capacity_NRHO) && ...
    all(load_L4   <= capacity_L4) && ...
    all(load_L5   <= capacity_L5);


%% ========================================================================
% Return
% ========================================================================

result.success = ...
    all(~isnan(selected_rank)) && ...
    capacity_feasible;

result.converged = ...
    converged;

result.rank = ...
    selected_rank;

result.GDOP = ...
    selected_GDOP;

result.NRHO_ID = ...
    selected_NRHO;

result.L4_ID = ...
    selected_L4;

result.L5_ID = ...
    selected_L5;

result.objective = ...
    sum(selected_GDOP);

result.load_NRHO = ...
    load_NRHO;

result.load_L4 = ...
    load_L4;

result.load_L5 = ...
    load_L5;

result.n_sweeps = ...
    actual_sweeps;

result.local_solve_count = ...
    local_solve_count;

result.objective_history = ...
    objective_history;

result.rank_history = ...
    rank_history;


if converged

    result.status = ...
        sprintf( ...
        'Converged after %d complete sweeps.', ...
        actual_sweeps);

else

    result.status = ...
        sprintf( ...
        'Maximum number of sweeps (%d) reached.', ...
        max_sweeps);

end

end


%% =========================================================================
% Local MILP
% =========================================================================

function [success,selected_rank,fval,exitflag] = ...
    solveLocalFRTMILP( ...
        i, ...
        G, ...
        NRHO, ...
        L4, ...
        L5, ...
        load_NRHO, ...
        load_L4, ...
        load_L5, ...
        capacity_NRHO, ...
        capacity_L4, ...
        capacity_L5, ...
        options)

%% ------------------------------------------------------------------------
% Local problem for FRT i:
%
%   min    sum_k G(i,k) x_i(k)
%
%   s.t.   sum_k x_i(k) = 1
%
%          resource usage <=
%          capacity - usage of all other FRTs
%
%          x_i(k) in {0,1}
%
% -------------------------------------------------------------------------

K = ...
    size(G,2);

N_NRHO = ...
    length(capacity_NRHO);

N_L4 = ...
    length(capacity_L4);

N_L5 = ...
    length(capacity_L5);


%% ------------------------------------------------------------------------
% Local binary variables
%
% x_i = [x_i1 ... x_iK]'
% -------------------------------------------------------------------------

nVar = ...
    K;

f = ...
    zeros(nVar,1);

lb = ...
    zeros(nVar,1);

ub = ...
    ones(nVar,1);


%% ------------------------------------------------------------------------
% Local resource constraint matrix
% -------------------------------------------------------------------------

A = ...
    zeros( ...
        N_NRHO + ...
        N_L4 + ...
        N_L5, ...
        K);


for k = 1:K

    n  = NRHO(i,k);
    l4 = L4(i,k);
    l5 = L5(i,k);


    valid_candidate = ...
        isfinite(G(i,k)) && ...
        isfinite(n) && ...
        isfinite(l4) && ...
        isfinite(l5) && ...
        n  >= 1 && n  <= N_NRHO && n  == round(n) && ...
        l4 >= 1 && l4 <= N_L4   && l4 == round(l4) && ...
        l5 >= 1 && l5 <= N_L5   && l5 == round(l5);


    if ~valid_candidate

        ub(k) = 0;
        continue;

    end


    %% Objective

    f(k) = ...
        G(i,k);


    %% NRHO incidence

    A(n,k) = ...
        1;


    %% L4 incidence

    A(N_NRHO + l4,k) = ...
        1;


    %% L5 incidence

    A(N_NRHO + N_L4 + l5,k) = ...
        1;

end


%% ------------------------------------------------------------------------
% Residual capacity
%
% load_* contains ONLY other FRT assignments.
% -------------------------------------------------------------------------

b = [
    capacity_NRHO(:) - load_NRHO(:)
    capacity_L4(:)   - load_L4(:)
    capacity_L5(:)   - load_L5(:)
    ];


%% ------------------------------------------------------------------------
% Exactly one candidate
% -------------------------------------------------------------------------

Aeq = ...
    ones(1,K);

beq = ...
    1;


%% ------------------------------------------------------------------------
% Binary variables
% -------------------------------------------------------------------------

intcon = ...
    1:K;


%% ------------------------------------------------------------------------
% Solve local MILP
% -------------------------------------------------------------------------

[x,fval,exitflag] = ...
    intlinprog( ...
        f, ...
        intcon, ...
        A,b, ...
        Aeq,beq, ...
        lb,ub, ...
        options);


%% ------------------------------------------------------------------------
% Decode local result
% -------------------------------------------------------------------------

success = ...
    exitflag > 0;

selected_rank = ...
    NaN;


if success

    selected_rank = ...
        find(x > 0.5,1);

    if isempty(selected_rank)

        success = false;
        selected_rank = NaN;

    end

end

end


%% =========================================================================
% Objective computation
% =========================================================================

function objective = ...
    computeAssignmentObjective( ...
        G, ...
        selected_rank)

N_FRT = ...
    size(G,1);

objective = ...
    0;


for i = 1:N_FRT

    objective = ...
        objective + ...
        G(i,selected_rank(i));

end

end