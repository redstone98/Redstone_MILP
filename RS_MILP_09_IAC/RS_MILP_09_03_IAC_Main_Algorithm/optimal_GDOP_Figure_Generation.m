function [t_vector] = optimal_GDOP_Figure_Generation(FRT_GDOP_lookup_table, t_vector)


%% ========================================================================
% Figure 1
% Physical FRT ID vs Time - Optimal GDOP Lookup Map
% ========================================================================

validIndex = ...
    isfinite(FRT_GDOP_lookup_table.BestGDOP);


figure( ...
    'Name','FRT ID-Aware GDOP Lookup Map');

scatter( ...
    FRT_GDOP_lookup_table.Time(validIndex), ...
    FRT_GDOP_lookup_table.FRT_ID(validIndex), ...
    35, ...
    FRT_GDOP_lookup_table.BestGDOP(validIndex), ...
    'filled');

xlabel('Time [TU]');
ylabel('Physical FRT Satellite ID');

title( ...
    'Optimal GDOP Lookup Map for Dynamic FRT Constellation');

cb = colorbar;
cb.Label.String = 'Minimum GDOP';

grid on;
box on;

xlim([t_vector(1) t_vector(end)]);

yticks( ...
    unique(FRT_GDOP_lookup_table.FRT_ID));


%% ========================================================================
% Figure 2
% Best Combination Index Lookup Map
% ========================================================================

figure( ...
    'Name','Best Combination Lookup Map');

scatter( ...
    FRT_GDOP_lookup_table.Time, ...
    FRT_GDOP_lookup_table.FRT_ID, ...
    35, ...
    FRT_GDOP_lookup_table.CombinationIndex, ...
    'filled');

xlabel('Time [TU]');
ylabel('Physical FRT Satellite ID');

title( ...
    'Best NRHO-L4-L5 Combination Lookup Map');

cb = colorbar;
cb.Label.String = 'Combination Index';

grid on;
box on;

xlim([t_vector(1) t_vector(end)]);

yticks( ...
    unique(FRT_GDOP_lookup_table.FRT_ID));


%% ========================================================================
% Figure 3
% Individual Servicing Satellite ID Lookup Maps
% ========================================================================

figure( ...
    'Name','Optimal Servicing Satellite IDs');

tiledlayout(3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
% NRHO
% -------------------------------------------------------------------------

nexttile;

scatter( ...
    FRT_GDOP_lookup_table.Time, ...
    FRT_GDOP_lookup_table.FRT_ID, ...
    30, ...
    FRT_GDOP_lookup_table.BestNRHO_ID, ...
    'filled');

ylabel('FRT ID');

title('Optimal NRHO Satellite');

cb = colorbar;
cb.Label.String = 'NRHO ID';

grid on;
box on;


%% ------------------------------------------------------------------------
% L4
% -------------------------------------------------------------------------

nexttile;

scatter( ...
    FRT_GDOP_lookup_table.Time, ...
    FRT_GDOP_lookup_table.FRT_ID, ...
    30, ...
    FRT_GDOP_lookup_table.BestL4_ID, ...
    'filled');

ylabel('FRT ID');

title('Optimal L4 Satellite');

cb = colorbar;
cb.Label.String = 'L4 ID';

grid on;
box on;


%% ------------------------------------------------------------------------
% L5
% -------------------------------------------------------------------------

nexttile;

scatter( ...
    FRT_GDOP_lookup_table.Time, ...
    FRT_GDOP_lookup_table.FRT_ID, ...
    30, ...
    FRT_GDOP_lookup_table.BestL5_ID, ...
    'filled');

xlabel('Time [TU]');
ylabel('FRT ID');

title('Optimal L5 Satellite');

cb = colorbar;
cb.Label.String = 'L5 ID';

grid on;
box on;


%% ========================================================================
% Figure 4
% GDOP History for Each Physical FRT Satellite
% ========================================================================

unique_FRT_IDs = ...
    unique(FRT_GDOP_lookup_table.FRT_ID);

figure( ...
    'Name','GDOP History by Physical FRT ID');

hold on;


for ii = 1:length(unique_FRT_IDs)

    current_ID = ...
        unique_FRT_IDs(ii);

    idx = ...
        FRT_GDOP_lookup_table.FRT_ID == current_ID;

    current_time = ...
        FRT_GDOP_lookup_table.Time(idx);

    current_GDOP = ...
        FRT_GDOP_lookup_table.BestGDOP(idx);

    [current_time,sort_index] = ...
        sort(current_time);

    current_GDOP = ...
        current_GDOP(sort_index);


    plot( ...
        current_time, ...
        current_GDOP, ...
        'LineWidth',1.4, ...
        'DisplayName', ...
        sprintf('FRT %d',current_ID));

end


xlabel('Time [TU]');
ylabel('Minimum GDOP');

title( ...
    'Optimal GDOP History of Individual FRT Satellites');

grid on;
box on;

legend( ...
    'Location','eastoutside');

hold off;



end