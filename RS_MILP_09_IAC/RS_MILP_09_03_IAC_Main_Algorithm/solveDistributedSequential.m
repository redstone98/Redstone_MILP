function result = solveDistributedSequential( ...
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
% Sequential Distributed Assignment
% ========================================================================

N_FRT = ...
    size(topK_GDOP_history,2);

K = ...
    size(topK_GDOP_history,3);


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


%% Resource load

load_NRHO = zeros(length(capacity_NRHO),1);
load_L4   = zeros(length(capacity_L4),1);
load_L5   = zeros(length(capacity_L5),1);


%% Results

selected_rank = NaN(N_FRT,1);

selected_GDOP = NaN(N_FRT,1);

selected_NRHO = NaN(N_FRT,1);
selected_L4   = NaN(N_FRT,1);
selected_L5   = NaN(N_FRT,1);


%% ========================================================================
% Sequential decision
% ========================================================================

for order_index = 1:N_FRT

    i = ...
        decision_order(order_index);


    for k = 1:K

        if ~isfinite(G(i,k))
            continue;
        end


        n = NRHO(i,k);
        l4 = L4(i,k);
        l5 = L5(i,k);


        feasible = ...
            load_NRHO(n) < capacity_NRHO(n) && ...
            load_L4(l4) < capacity_L4(l4) && ...
            load_L5(l5) < capacity_L5(l5);


        if feasible

            selected_rank(i) = k;

            selected_GDOP(i) = ...
                G(i,k);

            selected_NRHO(i) = n;
            selected_L4(i)   = l4;
            selected_L5(i)   = l5;


            load_NRHO(n) = ...
                load_NRHO(n) + 1;

            load_L4(l4) = ...
                load_L4(l4) + 1;

            load_L5(l5) = ...
                load_L5(l5) + 1;

            break;

        end

    end

end


%% ========================================================================
% Return
% ========================================================================

result.success = ...
    all(~isnan(selected_rank));

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
    sum(selected_GDOP,'omitnan');

result.load_NRHO = ...
    load_NRHO;

result.load_L4 = ...
    load_L4;

result.load_L5 = ...
    load_L5;

end