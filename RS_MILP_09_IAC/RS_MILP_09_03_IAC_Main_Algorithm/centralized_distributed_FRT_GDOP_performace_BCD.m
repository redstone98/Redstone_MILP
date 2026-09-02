function [central_assignment_table,  distributed_assignment_table, FRT_GDOP_performance_table, ...
 central_success_history, central_objective_history , central_runtime_history, central_rank_history, central_GDOP_history, central_NRHO_history, central_L4_history, central_L5_history, ...
 distributed_success_history, distributed_objective_history , distributed_runtime_history, distributed_rank_history, distributed_GDOP_history, distributed_NRHO_history, distributed_L4_history, distributed_L5_history] = ...
 centralized_distributed_FRT_GDOP_performace_BCD(t_vector, topK_GDOP_history, ...
                                                            topK_NRHO_ID_history, ...
                                                            topK_L4_ID_history, ...
                                                            topK_L5_ID_history, ...
                                                            capacity_NRHO, ...
                                                            capacity_L4, ...
                                                            capacity_L5, ...
                                                            FRT_ID_history, ...
                                                            FRT_age_history)


%% ========================================================================
% VI. ALL-TIME CENTRALIZED / DISTRIBUTED ASSIGNMENT
% ========================================================================

nTime = size(topK_GDOP_history,1);
N_FRT = size(topK_GDOP_history,2);

%% ------------------------------------------------------------------------
% Capacity
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% Preallocation - Centralized
% -------------------------------------------------------------------------

central_success_history   = false(nTime,1);
central_objective_history = NaN(nTime,1);
central_runtime_history   = NaN(nTime,1);

central_rank_history      = NaN(nTime,N_FRT);
central_GDOP_history      = NaN(nTime,N_FRT);

central_NRHO_history      = NaN(nTime,N_FRT);
central_L4_history        = NaN(nTime,N_FRT);
central_L5_history        = NaN(nTime,N_FRT);


%% ------------------------------------------------------------------------
% Preallocation - Distributed Sequential
% -------------------------------------------------------------------------

distributed_success_history   = false(nTime,1);
distributed_objective_history = NaN(nTime,1);
distributed_runtime_history   = NaN(nTime,1);

distributed_rank_history      = NaN(nTime,N_FRT);
distributed_GDOP_history      = NaN(nTime,N_FRT);

distributed_NRHO_history      = NaN(nTime,N_FRT);
distributed_L4_history        = NaN(nTime,N_FRT);
distributed_L5_history        = NaN(nTime,N_FRT);


%% ========================================================================
% Solve all timesteps
% ========================================================================

for time_index = 1:nTime

    % fprintf( ...
    %     '\n====================================================\n');
    % fprintf( ...
    %     'Assignment timestep %d / %d\n', ...
    %     time_index, nTime);
    % fprintf( ...
    %     '====================================================\n');


    %% ====================================================================
    % 1. Centralized
    % ====================================================================

    tic;

    central_i = solveCentralizedGDOPAssignment( ...
        time_index, ...
        topK_GDOP_history, ...
        topK_NRHO_ID_history, ...
        topK_L4_ID_history, ...
        topK_L5_ID_history, ...
        capacity_NRHO, ...
        capacity_L4, ...
        capacity_L5);

    central_runtime_history(time_index) = toc;


    %% Store success / objective

    central_success_history(time_index) = ...
        central_i.success;

    central_objective_history(time_index) = ...
        central_i.objective;


    %% Store satellite-level results

    if central_i.success

        central_rank_history(time_index,:) = ...
            central_i.rank';

        central_GDOP_history(time_index,:) = ...
            central_i.GDOP';

        central_NRHO_history(time_index,:) = ...
            central_i.NRHO_ID';

        central_L4_history(time_index,:) = ...
            central_i.L4_ID';

        central_L5_history(time_index,:) = ...
            central_i.L5_ID';

    end


    %% ====================================================================
    % 2. Distributed Sequential
    % ====================================================================

    decision_order = ...
        1:N_FRT;

    tic;

    distributed_i = solveDistributedSequential_BCD( ...
        time_index, ...
        topK_GDOP_history, ...
        topK_NRHO_ID_history, ...
        topK_L4_ID_history, ...
        topK_L5_ID_history, ...
        capacity_NRHO, ...
        capacity_L4, ...
        capacity_L5, ...
        decision_order);

    distributed_runtime_history(time_index) = toc;


    %% Store success / objective

    distributed_success_history(time_index) = ...
        distributed_i.success;

    distributed_objective_history(time_index) = ...
        distributed_i.objective;


    %% Store satellite-level results

    if distributed_i.success

        distributed_rank_history(time_index,:) = ...
            distributed_i.rank';

        distributed_GDOP_history(time_index,:) = ...
            distributed_i.GDOP';

        distributed_NRHO_history(time_index,:) = ...
            distributed_i.NRHO_ID';

        distributed_L4_history(time_index,:) = ...
            distributed_i.L4_ID';

        distributed_L5_history(time_index,:) = ...
            distributed_i.L5_ID';

    end

end

%% ========================================================================
% VII. SUCCESS CHECK
% ========================================================================

all_central_success = ...
    all(central_success_history);

all_distributed_success = ...
    all(distributed_success_history);

all_methods_success = ...
    all_central_success && ...
    all_distributed_success;


fprintf('\n');
fprintf('====================================================\n');
fprintf('ALL-TIME ASSIGNMENT SUCCESS SUMMARY\n');
fprintf('====================================================\n');

fprintf( ...
    'Centralized Success   : %d / %d\n', ...
    sum(central_success_history), ...
    nTime);

fprintf( ...
    'Distributed Success   : %d / %d\n', ...
    sum(distributed_success_history), ...
    nTime);


if all_central_success
    fprintf('Centralized : ALL TIMESTEPS SUCCESS\n');
else
    fprintf('Centralized : FAILURE EXISTS\n');
end


if all_distributed_success
    fprintf('Distributed : ALL TIMESTEPS SUCCESS\n');
else
    fprintf('Distributed : FAILURE EXISTS\n');
end

%% ========================================================================
% VIII. STACK ALL ASSIGNMENT RESULTS
% ========================================================================

%% ------------------------------------------------------------------------
% Common indexing vectors
% -------------------------------------------------------------------------

TimeIndex = ...
    repmat((1:nTime)',N_FRT,1);

Time = ...
    repmat(t_vector(:),N_FRT,1);

FRTSlot = ...
    repelem((1:N_FRT)',nTime);

FRT_ID = ...
    reshape(FRT_ID_history,[],1);

FRT_PhaseAge = ...
    reshape(FRT_age_history,[],1);


%% ========================================================================
% Centralized Table
% ========================================================================

Central_Success = ...
    repmat(central_success_history,N_FRT,1);

Central_Rank = ...
    reshape(central_rank_history,[],1);

Central_GDOP = ...
    reshape(central_GDOP_history,[],1);

Central_NRHO_ID = ...
    reshape(central_NRHO_history,[],1);

Central_L4_ID = ...
    reshape(central_L4_history,[],1);

Central_L5_ID = ...
    reshape(central_L5_history,[],1);


central_assignment_table = table( ...
    TimeIndex, ...
    Time, ...
    FRTSlot, ...
    FRT_ID, ...
    FRT_PhaseAge, ...
    Central_Success, ...
    Central_Rank, ...
    Central_GDOP, ...
    Central_NRHO_ID, ...
    Central_L4_ID, ...
    Central_L5_ID);


%% ========================================================================
% Distributed Table
% ========================================================================

Distributed_Success = ...
    repmat(distributed_success_history,N_FRT,1);

Distributed_Rank = ...
    reshape(distributed_rank_history,[],1);

Distributed_GDOP = ...
    reshape(distributed_GDOP_history,[],1);

Distributed_NRHO_ID = ...
    reshape(distributed_NRHO_history,[],1);

Distributed_L4_ID = ...
    reshape(distributed_L4_history,[],1);

Distributed_L5_ID = ...
    reshape(distributed_L5_history,[],1);


distributed_assignment_table = table( ...
    TimeIndex, ...
    Time, ...
    FRTSlot, ...
    FRT_ID, ...
    FRT_PhaseAge, ...
    Distributed_Success, ...
    Distributed_Rank, ...
    Distributed_GDOP, ...
    Distributed_NRHO_ID, ...
    Distributed_L4_ID, ...
    Distributed_L5_ID);

%% ========================================================================
% IX. GDOP PERFORMANCE BY PHYSICAL FRT SATELLITE
% ========================================================================

if all_methods_success

    fprintf('\n');
    fprintf('All timesteps were successfully assigned.\n');
    fprintf('Calculating GDOP performance by physical FRT ID...\n');


    %% --------------------------------------------------------------------
    % All physical FRT IDs appearing during the simulation
    % ---------------------------------------------------------------------

    unique_FRT_IDs = ...
        unique(FRT_ID_history(:));

    nPhysicalFRT = ...
        length(unique_FRT_IDs);


    %% --------------------------------------------------------------------
    % Preallocation
    % ---------------------------------------------------------------------

    ServiceCount = zeros(nPhysicalFRT,1);


    %% Centralized statistics

    Central_MeanGDOP   = NaN(nPhysicalFRT,1);
    Central_MedianGDOP = NaN(nPhysicalFRT,1);
    Central_MinGDOP    = NaN(nPhysicalFRT,1);
    Central_MaxGDOP    = NaN(nPhysicalFRT,1);
    Central_StdGDOP    = NaN(nPhysicalFRT,1);


    %% Distributed statistics

    Distributed_MeanGDOP   = NaN(nPhysicalFRT,1);
    Distributed_MedianGDOP = NaN(nPhysicalFRT,1);
    Distributed_MinGDOP    = NaN(nPhysicalFRT,1);
    Distributed_MaxGDOP    = NaN(nPhysicalFRT,1);
    Distributed_StdGDOP    = NaN(nPhysicalFRT,1);


    %% ====================================================================
    % Physical ID loop
    % ====================================================================

    for ii = 1:nPhysicalFRT

        current_ID = ...
            unique_FRT_IDs(ii);


        %% -------------------------------------------------------------
        % Centralized
        % -------------------------------------------------------------

        idx_central = ...
            central_assignment_table.FRT_ID == current_ID;

        gdop_central = ...
            central_assignment_table.Central_GDOP(idx_central);

        gdop_central = ...
            gdop_central(isfinite(gdop_central));


        %% -------------------------------------------------------------
        % Distributed
        % -------------------------------------------------------------

        idx_distributed = ...
            distributed_assignment_table.FRT_ID == current_ID;

        gdop_distributed = ...
            distributed_assignment_table.Distributed_GDOP( ...
                idx_distributed);

        gdop_distributed = ...
            gdop_distributed(isfinite(gdop_distributed));


        %% -------------------------------------------------------------
        % Service count
        % -------------------------------------------------------------

        ServiceCount(ii) = ...
            length(gdop_central);


        %% -------------------------------------------------------------
        % Centralized statistics
        % -------------------------------------------------------------

        if ~isempty(gdop_central)

            Central_MeanGDOP(ii) = ...
                mean(gdop_central);

            Central_MedianGDOP(ii) = ...
                median(gdop_central);

            Central_MinGDOP(ii) = ...
                min(gdop_central);

            Central_MaxGDOP(ii) = ...
                max(gdop_central);

            Central_StdGDOP(ii) = ...
                std(gdop_central);

        end


        %% -------------------------------------------------------------
        % Distributed statistics
        % -------------------------------------------------------------

        if ~isempty(gdop_distributed)

            Distributed_MeanGDOP(ii) = ...
                mean(gdop_distributed);

            Distributed_MedianGDOP(ii) = ...
                median(gdop_distributed);

            Distributed_MinGDOP(ii) = ...
                min(gdop_distributed);

            Distributed_MaxGDOP(ii) = ...
                max(gdop_distributed);

            Distributed_StdGDOP(ii) = ...
                std(gdop_distributed);

        end

    end


    %% ====================================================================
    % Comparison
    % ====================================================================

    MeanGDOP_Difference = ...
        Distributed_MeanGDOP - ...
        Central_MeanGDOP;

    MeanGDOP_Degradation_percent = ...
        100 * ...
        MeanGDOP_Difference ./ ...
        Central_MeanGDOP;


    %% ====================================================================
    % Final performance table
    % ====================================================================

    FRT_GDOP_performance_table = table( ...
        unique_FRT_IDs, ...
        ServiceCount, ...
        Central_MeanGDOP, ...
        Central_MedianGDOP, ...
        Central_MinGDOP, ...
        Central_MaxGDOP, ...
        Central_StdGDOP, ...
        Distributed_MeanGDOP, ...
        Distributed_MedianGDOP, ...
        Distributed_MinGDOP, ...
        Distributed_MaxGDOP, ...
        Distributed_StdGDOP, ...
        MeanGDOP_Difference, ...
        MeanGDOP_Degradation_percent, ...
        'VariableNames',{ ...
            'FRT_ID', ...
            'ServiceEpochCount', ...
            'Central_MeanGDOP', ...
            'Central_MedianGDOP', ...
            'Central_MinGDOP', ...
            'Central_MaxGDOP', ...
            'Central_StdGDOP', ...
            'Distributed_MeanGDOP', ...
            'Distributed_MedianGDOP', ...
            'Distributed_MinGDOP', ...
            'Distributed_MaxGDOP', ...
            'Distributed_StdGDOP', ...
            'MeanGDOP_Difference', ...
            'MeanGDOP_Degradation_Percent'});


    disp(FRT_GDOP_performance_table);


else

    warning( ...
        ['At least one assignment method failed at one or more ', ...
         'timesteps. GDOP performance table was not generated.']);

end


end