function [min_GDOP, best_combination_index, best_NRHO_index, best_L4_index, best_L5_index] = ...
    get_best_GDOP_info(rv_FRT_one, GDOP_history, NRHO_index_history, L4_index_history, L5_index_history, t_vector, number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs)


%% ========================================================================
%  FRT Position-Based Best Combination Lookup Chart
% ========================================================================

% Best combination at each timestep
[min_GDOP, best_combination_index] = min(GDOP_history, [], 2);

best_NRHO_index = NRHO_index_history(best_combination_index);
best_L4_index   = L4_index_history(best_combination_index);
best_L5_index   = L5_index_history(best_combination_index);

%% ------------------------------------------------------------------------
% FRT trajectory information
% -------------------------------------------------------------------------

r_FRT_history = rv_FRT_one(:,1:3);

% Distance traveled along the FRT trajectory
FRT_path_distance = zeros(length(t_vector),1);

for time_index = 2:length(t_vector)

    dr = ...
        r_FRT_history(time_index,:) - ...
        r_FRT_history(time_index-1,:);

    FRT_path_distance(time_index) = ...
        FRT_path_distance(time_index-1) + norm(dr);

end

% Normalized orbit position
% 0 = initial FRT position
% 1 = final FRT position
FRT_orbit_position = ...
    FRT_path_distance ./ FRT_path_distance(end);


%% ========================================================================
% Figure 1
% Best NRHO-L4-L5 Combination Lookup Chart
% ========================================================================

figure('Name','Best Navigation Satellite Combination Lookup Chart');

tiledlayout(4,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
% 1. Best NRHO satellite index
% -------------------------------------------------------------------------

nexttile;

stairs( ...
    FRT_orbit_position, ...
    best_NRHO_index, ...
    'LineWidth',1.5);

ylabel('NRHO Index');

title('Best NRHO-L4-L5 Combination Along FRT Trajectory');

grid on;
box on;

ylim([0.5 number_of_NRHO_SATs + 0.5]);
yticks(1:number_of_NRHO_SATs);


%% ------------------------------------------------------------------------
% 2. Best L4 satellite index
% -------------------------------------------------------------------------

nexttile;

stairs( ...
    FRT_orbit_position, ...
    best_L4_index, ...
    'LineWidth',1.5);

ylabel('L4 Index');

grid on;
box on;

ylim([0.5 number_of_L4_SATs + 0.5]);
yticks(1:number_of_L4_SATs);


%% ------------------------------------------------------------------------
% 3. Best L5 satellite index
% -------------------------------------------------------------------------

nexttile;

stairs( ...
    FRT_orbit_position, ...
    best_L5_index, ...
    'LineWidth',1.5);

ylabel('L5 Index');

grid on;
box on;

ylim([0.5 number_of_L5_SATs + 0.5]);
yticks(1:number_of_L5_SATs);


%% ------------------------------------------------------------------------
% 4. Best combination index
% -------------------------------------------------------------------------

nexttile;

stairs( ...
    FRT_orbit_position, ...
    best_combination_index, ...
    'LineWidth',1.5);

xlabel('Normalized FRT Orbit Position');
ylabel('Combination Index');

grid on;
box on;

xlim([0 1]);

sgtitle( ...
    'FRT Orbit Position-Based Navigation Satellite Lookup Chart');


%% ========================================================================
% Figure 2
% Minimum GDOP vs Time
% ========================================================================

figure('Name','Minimum GDOP History');

plot( ...
    t_vector, ...
    min_GDOP, ...
    'LineWidth',1.8);

xlabel('Time [TU]');
ylabel('Minimum GDOP');

title('Minimum Achievable GDOP Along FRT Trajectory');

grid on;
box on;

xlim([t_vector(1) t_vector(end)]);

end