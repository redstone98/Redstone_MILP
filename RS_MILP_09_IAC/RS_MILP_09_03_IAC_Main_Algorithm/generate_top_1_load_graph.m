function [max_NRHO_load, max_L4_load, max_L5_load] = ...
    generate_top_1_load_graph(number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs, ...
                              NRHO_load_history, L4_load_history, L5_load_history, t_vector)


N_NRHO = number_of_NRHO_SATs;
N_L4   = number_of_L4_SATs;
N_L5   = number_of_L5_SATs;

%% ========================================================================
% Maximum servicing burden
% ========================================================================

[max_NRHO_load, linear_idx] = ...
    max(NRHO_load_history(:));

[max_NRHO_time, max_NRHO_ID] = ...
    ind2sub(size(NRHO_load_history),linear_idx);


[max_L4_load, linear_idx] = ...
    max(L4_load_history(:));

[max_L4_time, max_L4_ID] = ...
    ind2sub(size(L4_load_history),linear_idx);


[max_L5_load, linear_idx] = ...
    max(L5_load_history(:));

[max_L5_time, max_L5_ID] = ...
    ind2sub(size(L5_load_history),linear_idx);


fprintf('\n========================================\n');
fprintf('Maximum Top-1 Servicing Load\n');
fprintf('========================================\n');

fprintf( ...
    'NRHO %d : %d FRT satellites at timestep %d\n', ...
    max_NRHO_ID, ...
    max_NRHO_load, ...
    max_NRHO_time);

fprintf( ...
    'L4   %d : %d FRT satellites at timestep %d\n', ...
    max_L4_ID, ...
    max_L4_load, ...
    max_L4_time);

fprintf( ...
    'L5   %d : %d FRT satellites at timestep %d\n', ...
    max_L5_ID, ...
    max_L5_load, ...
    max_L5_time);

%% ========================================================================
% Top-1 Load Heatmap
% ========================================================================

figure( ...
    'Name','Top-1 Servicing Satellite Load');

tiledlayout(3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% NRHO

nexttile;

imagesc( ...
    t_vector, ...
    1:N_NRHO, ...
    NRHO_load_history');

axis xy;

ylabel('NRHO Satellite ID');

title( ...
    'NRHO Top-1 Request Load');

cb = colorbar;
cb.Label.String = 'Number of FRT Requests';


%% L4

nexttile;

imagesc( ...
    t_vector, ...
    1:N_L4, ...
    L4_load_history');

axis xy;

ylabel('L4 Satellite ID');

title( ...
    'L4 Top-1 Request Load');

cb = colorbar;
cb.Label.String = 'Number of FRT Requests';


%% L5

nexttile;

imagesc( ...
    t_vector, ...
    1:N_L5, ...
    L5_load_history');

axis xy;

xlabel('Time [TU]');
ylabel('L5 Satellite ID');

title( ...
    'L5 Top-1 Request Load');

cb = colorbar;
cb.Label.String = 'Number of FRT Requests';

end