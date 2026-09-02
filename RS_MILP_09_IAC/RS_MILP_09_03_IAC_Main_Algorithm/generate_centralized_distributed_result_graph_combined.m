function [unique_FRT_IDs] = generate_centralized_distributed_result_graph(central_assignment_table,  distributed_assignment_table, t_vector)


%% ========================================================================
% Individual Physical FRT Satellite GDOP History
% CENTRALIZED ASSIGNMENT
% ========================================================================

unique_FRT_IDs = ...
    unique(central_assignment_table.FRT_ID);

figure( ...
    'Name','Centralized, Distributed GDOP History - Individual FRT Satellites');

hold on;


for ii = 1:length(unique_FRT_IDs)

    current_ID = ...
        unique_FRT_IDs(ii);

    %% Extract current physical FRT satellite
    idx = ...
        central_assignment_table.FRT_ID == current_ID & ...
        central_assignment_table.Central_Success;

    current_time = ...
        central_assignment_table.Time(idx);

    current_GDOP = ...
        central_assignment_table.Central_GDOP(idx);


    %% Remove invalid points
    valid_idx = ...
        isfinite(current_GDOP);

    current_time = ...
        current_time(valid_idx);

    current_GDOP = ...
        current_GDOP(valid_idx);


    %% Sort by time
    [current_time,sort_idx] = ...
        sort(current_time);

    current_GDOP = ...
        current_GDOP(sort_idx);


    %% Plot
    plot( ...
        current_time, ...
        current_GDOP, ...
        'LineWidth',1, ...
        'Color','blue');

end


xlabel('Time [TU]');
ylabel('Assigned GDOP');

title( ...
    ' GDOP History of Individual FRT Satellites, Blue: Centralized, Red: Distributed');

xlim([t_vector(1) t_vector(end)]);

grid on;
box on;




%% ========================================================================
% Individual Physical FRT Satellite GDOP History
% DISTRIBUTED SEQUENTIAL ASSIGNMENT
% ========================================================================

unique_FRT_IDs = ...
    unique(distributed_assignment_table.FRT_ID);



for ii = 1:length(unique_FRT_IDs)

    current_ID = ...
        unique_FRT_IDs(ii);

    %% Extract current physical FRT satellite
    idx = ...
        distributed_assignment_table.FRT_ID == current_ID & ...
        distributed_assignment_table.Distributed_Success;

    current_time = ...
        distributed_assignment_table.Time(idx);

    current_GDOP = ...
        distributed_assignment_table.Distributed_GDOP(idx);


    %% Remove invalid points
    valid_idx = ...
        isfinite(current_GDOP);

    current_time = ...
        current_time(valid_idx);

    current_GDOP = ...
        current_GDOP(valid_idx);


    %% Sort by time
    [current_time,sort_idx] = ...
        sort(current_time);

    current_GDOP = ...
        current_GDOP(sort_idx);


    %% Plot
    plot( ...
        current_time, ...
        current_GDOP, ...
        'LineWidth',1.5, ...
        'Color','r', 'LineStyle','--');

end


hold off;

end