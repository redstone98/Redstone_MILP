function [central_GDOP_plot] = plot_centralized_distributed_graph(t_vector, central_assignment_table,  distributed_assignment_table, FRT_GDOP_performance_table, ...
 central_success_history, central_objective_history , central_runtime_history, central_rank_history, ...
 central_GDOP_history, central_NRHO_history, central_L4_history, central_L5_history, ...
 distributed_success_history, distributed_objective_history , distributed_runtime_history, distributed_rank_history, ...
 distributed_GDOP_history, distributed_NRHO_history, distributed_L4_history, distributed_L5_history, ...
 topK_GDOP_history)


%% ========================================================================
% XIV. PLOT ALL RESULTS
% ========================================================================

%% ------------------------------------------------------------------------
% Remove Inf for plotting/statistics
% -------------------------------------------------------------------------

central_GDOP_plot = central_GDOP_history;
distributed_GDOP_plot = distributed_GDOP_history;

central_GDOP_plot(~isfinite(central_GDOP_plot)) = NaN;
distributed_GDOP_plot(~isfinite(distributed_GDOP_plot)) = NaN;


%% ------------------------------------------------------------------------
% Time-wise statistics
% -------------------------------------------------------------------------

central_mean_GDOP_time = ...
    mean(central_GDOP_plot,2,'omitnan');

central_max_GDOP_time = ...
    max(central_GDOP_plot,[],2,'omitnan');

central_min_GDOP_time = ...
    min(central_GDOP_plot,[],2,'omitnan');


distributed_mean_GDOP_time = ...
    mean(distributed_GDOP_plot,2,'omitnan');

distributed_max_GDOP_time = ...
    max(distributed_GDOP_plot,[],2,'omitnan');

distributed_min_GDOP_time = ...
    min(distributed_GDOP_plot,[],2,'omitnan');


%% ========================================================================
% FIGURE 1
% Assignment Success History
% ========================================================================

figure( ...
    'Name','Assignment Success History');

hold on;

stairs( ...
    t_vector, ...
    double(central_success_history), ...
    'LineWidth',1.8);

stairs( ...
    t_vector, ...
    double(distributed_success_history), ...
    '--', ...
    'LineWidth',1.8);

xlabel('Time [TU]');
ylabel('Assignment Success');

title('Assignment Success History');

legend( ...
    'Centralized', ...
    'Distributed Sequential', ...
    'Location','best');

yticks([0 1]);
yticklabels({'Failure','Success'});

ylim([-0.1 1.1]);

grid on;
box on;

hold off;


%% ========================================================================
% FIGURE 2
% Mean / Worst-Case GDOP vs Time
% ========================================================================

figure( ...
    'Name','GDOP Performance vs Time');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
% Mean GDOP
% -------------------------------------------------------------------------

nexttile;

hold on;

plot( ...
    t_vector, ...
    central_mean_GDOP_time, ...
    'LineWidth',1.8);

plot( ...
    t_vector, ...
    distributed_mean_GDOP_time, ...
    '--', ...
    'LineWidth',1.8);

ylabel('Mean GDOP');

title('Mean GDOP of Active FRT Satellites');

legend( ...
    'Centralized', ...
    'Distributed Sequential', ...
    'Location','best');

grid on;
box on;

hold off;


%% ------------------------------------------------------------------------
% Worst-case GDOP
% -------------------------------------------------------------------------

nexttile;

hold on;

plot( ...
    t_vector, ...
    central_max_GDOP_time, ...
    'LineWidth',1.8);

plot( ...
    t_vector, ...
    distributed_max_GDOP_time, ...
    '--', ...
    'LineWidth',1.8);

xlabel('Time [TU]');
ylabel('Worst-Case GDOP');

title('Worst-Case GDOP Among Active FRT Satellites');

legend( ...
    'Centralized', ...
    'Distributed Sequential', ...
    'Location','best');

grid on;
box on;

hold off;


%% ========================================================================
% FIGURE 3
% Physical FRT Satellite ID - GDOP Lookup Map
% ========================================================================

all_methods_success = true;

if all_methods_success

    figure( ...
        'Name','Physical FRT GDOP Lookup Map');

    tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    %% --------------------------------------------------------------------
    % Centralized
    % ---------------------------------------------------------------------

    nexttile;

    valid_idx = ...
        isfinite(central_assignment_table.Central_GDOP);

    scatter( ...
        central_assignment_table.Time(valid_idx), ...
        central_assignment_table.FRT_ID(valid_idx), ...
        40, ...
        central_assignment_table.Central_GDOP(valid_idx), ...
        'filled');

    ylabel('Physical FRT Satellite ID');

    title('Centralized Assignment');

    cb = colorbar;
    cb.Label.String = 'GDOP';

    grid on;
    box on;


    %% --------------------------------------------------------------------
    % Distributed
    % ---------------------------------------------------------------------

    nexttile;

    valid_idx = ...
        isfinite(distributed_assignment_table.Distributed_GDOP);

    scatter( ...
        distributed_assignment_table.Time(valid_idx), ...
        distributed_assignment_table.FRT_ID(valid_idx), ...
        40, ...
        distributed_assignment_table.Distributed_GDOP(valid_idx), ...
        'filled');

    xlabel('Time [TU]');
    ylabel('Physical FRT Satellite ID');

    title('Distributed Sequential Assignment');

    cb = colorbar;
    cb.Label.String = 'GDOP';

    grid on;
    box on;

end


%% ========================================================================
% FIGURE 4
% GDOP Performance by Physical FRT Satellite
% ========================================================================



if all_methods_success

    FRT_ID_plot = ...
        FRT_GDOP_performance_table.FRT_ID;


    figure( ...
        'Name','GDOP by Physical FRT Satellite');

    tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    %% --------------------------------------------------------------------
    % Mean GDOP
    % ---------------------------------------------------------------------

    nexttile;

    Y_mean = [
        FRT_GDOP_performance_table.Central_MeanGDOP, ...
        FRT_GDOP_performance_table.Distributed_MeanGDOP
        ];

    bar( ...
        FRT_ID_plot, ...
        Y_mean);

    ylabel('Mean GDOP');

    title('Mean GDOP by Physical FRT Satellite');

    legend( ...
        'Centralized', ...
        'Distributed Sequential', ...
        'Location','best');

    xticks(FRT_ID_plot);

    grid on;
    box on;


    %% --------------------------------------------------------------------
    % Maximum GDOP
    % ---------------------------------------------------------------------

    nexttile;

    Y_max = [
        FRT_GDOP_performance_table.Central_MaxGDOP, ...
        FRT_GDOP_performance_table.Distributed_MaxGDOP
        ];

    bar( ...
        FRT_ID_plot, ...
        Y_max);

    xlabel('Physical FRT Satellite ID');
    ylabel('Maximum GDOP');

    title('Maximum GDOP by Physical FRT Satellite');

    legend( ...
        'Centralized', ...
        'Distributed Sequential', ...
        'Location','best');

    xticks(FRT_ID_plot);

    grid on;
    box on;

end


%% ========================================================================
% FIGURE 5
% Distributed GDOP Degradation by Physical FRT Satellite
% ========================================================================

if all_methods_success

    figure( ...
        'Name','Distributed GDOP Degradation');

    bar( ...
        FRT_GDOP_performance_table.FRT_ID, ...
        FRT_GDOP_performance_table.MeanGDOP_Degradation_Percent);

    xlabel('Physical FRT Satellite ID');
    ylabel('Mean GDOP Degradation [%]');

    title( ...
        'Distributed Sequential GDOP Degradation Relative to Centralized');

    yline( ...
        0, ...
        '--');

    xticks( ...
        FRT_GDOP_performance_table.FRT_ID);

    grid on;
    box on;

end


%% ========================================================================
% FIGURE 6
% Selected Candidate Rank Performance
% ========================================================================

if all_methods_success

    %% Number of Top-K candidates
    K = size(topK_GDOP_history,3);

    %% --------------------------------------------------------------------
    % Mean rank at each timestep
    % ---------------------------------------------------------------------

    central_mean_rank_time = ...
        mean(central_rank_history,2,'omitnan');

    distributed_mean_rank_time = ...
        mean(distributed_rank_history,2,'omitnan');


    figure( ...
        'Name','Candidate Rank Performance');

    tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    %% --------------------------------------------------------------------
    % Mean selected rank
    % ---------------------------------------------------------------------

    nexttile;

    hold on;

    plot( ...
        t_vector, ...
        central_mean_rank_time, ...
        'LineWidth',1.8);

    plot( ...
        t_vector, ...
        distributed_mean_rank_time, ...
        '--', ...
        'LineWidth',1.8);

    ylabel('Mean Selected Rank');

    title('Mean Top-K Candidate Rank');

    legend( ...
        'Centralized', ...
        'Distributed Sequential', ...
        'Location','best');

    grid on;
    box on;

    hold off;


    %% --------------------------------------------------------------------
    % Rank histogram
    % ---------------------------------------------------------------------

    nexttile;

    rank_axis = 1:K;

    central_rank_count = ...
        zeros(K,1);

    distributed_rank_count = ...
        zeros(K,1);


    for k = 1:K

        central_rank_count(k) = ...
            sum(central_rank_history(:) == k);

        distributed_rank_count(k) = ...
            sum(distributed_rank_history(:) == k);

    end


    %% --------------------------------------------------------------------
    % Percentage
    % ---------------------------------------------------------------------

    central_total = ...
        sum(central_rank_count);

    distributed_total = ...
        sum(distributed_rank_count);


    central_rank_percentage = ...
        100 * central_rank_count / central_total;

    distributed_rank_percentage = ...
        100 * distributed_rank_count / distributed_total;


    %% --------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------

    bar( ...
        rank_axis, ...
        [ ...
        central_rank_percentage, ...
        distributed_rank_percentage ...
        ]);

    xlabel('Selected Candidate Rank');
    ylabel('Selection Ratio [%]');

    title('Top-K Candidate Selection Distribution');

    legend( ...
        'Centralized', ...
        'Distributed Sequential', ...
        'Location','best');

    xticks(rank_axis);

    xlim([0.5 K+0.5]);

    grid on;
    box on;

end


%% ========================================================================
% FIGURE 7
% Runtime Comparison
% ========================================================================

figure( ...
    'Name','Computation Time Comparison');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
% Runtime over time
% -------------------------------------------------------------------------

nexttile;

hold on;

plot( ...
    t_vector, ...
    central_runtime_history, ...
    'LineWidth',1.6);

plot( ...
    t_vector, ...
    distributed_runtime_history, ...
    '--', ...
    'LineWidth',1.6);

ylabel('Runtime [sec]');

title('Assignment Computation Time');

legend( ...
    'Centralized', ...
    'Distributed Sequential', ...
    'Location','best');

grid on;
box on;

hold off;


%% ------------------------------------------------------------------------
% Mean runtime
% -------------------------------------------------------------------------

nexttile;

mean_runtime = [
    mean(central_runtime_history,'omitnan')
    mean(distributed_runtime_history,'omitnan')
    ];

bar( ...
    categorical({'Centralized','Distributed'}), ...
    mean_runtime);

ylabel('Mean Runtime [sec]');

title('Average Computation Time per Timestep');

grid on;
box on;


%% ========================================================================
% FIGURE 8
% Objective and Optimality Gap
% ========================================================================

if all_methods_success

    objective_gap_percent = ...
        100 * ...
        (distributed_objective_history - ...
         central_objective_history) ./ ...
        central_objective_history;


    figure( ...
        'Name','Objective Performance Comparison');

    tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    %% --------------------------------------------------------------------
    % Objective
    % ---------------------------------------------------------------------

    nexttile;

    hold on;

    plot( ...
        t_vector, ...
        central_objective_history, ...
        'LineWidth',1.8);

    plot( ...
        t_vector, ...
        distributed_objective_history, ...
        '--', ...
        'LineWidth',1.8);

    ylabel('\Sigma GDOP');

    title('Total GDOP Objective');

    legend( ...
        'Centralized', ...
        'Distributed Sequential', ...
        'Location','best');

    grid on;
    box on;

    hold off;


    %% --------------------------------------------------------------------
    % Optimality gap
    % ---------------------------------------------------------------------

    nexttile;

    plot( ...
        t_vector, ...
        objective_gap_percent, ...
        'LineWidth',1.8);

    yline( ...
        0, ...
        '--');

    xlabel('Time [TU]');
    ylabel('Optimality Gap [%]');

    title( ...
        'Distributed Optimality Gap Relative to Centralized');

    grid on;
    box on;

end




end
