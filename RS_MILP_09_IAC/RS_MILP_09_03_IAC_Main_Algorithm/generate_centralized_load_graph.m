function [max_NRHO_load, max_L4_load, max_L5_load] = ...
    generate_centralized_load_graph(number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs, ...
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
fprintf('Maximum Centralized Servicing Load\n');
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
% Centralized Load Heatmap
% ========================================================================

figure( ...
    'Name','Centralized Servicing Satellite Load');

tiledlayout(3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');
cmap = [...
    0.10 0.10 0.10;   % 0 : 거의 검정
    0.10 0.45 0.95;   % 1 : 선명한 파랑
    1.00 0.85 0.00;   % 2 : 노랑
    1.00 0.55 0.00;   % 3 : 주황
    0.35 0.75 0.25;   % 4 : 연두
    0.00 0.55 0.20;   % 5 : 진한 초록
    0.85 0.20 0.65;   % 6 : 자홍/마젠타
    0.85 0.00 0.00];  % 7 : 빨강


%% NRHO

nexttile;

imagesc( ...
    t_vector, ...
    1:N_NRHO, ...
    NRHO_load_history');

axis xy;

ylabel('NRHO Satellite ID');

title( ...
    'NRHO Centralized Request Load');

colormap(cmap);
clim([-0.5 7.5]);   % 0~7 정수값을 정확히 색 하나씩 대응
cb = colorbar;
cb.Ticks = 0:7;
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
    'L4 Centralized Request Load');

colormap(cmap);
clim([-0.5 7.5]);   % 0~7 정수값을 정확히 색 하나씩 대응
cb = colorbar;
cb.Ticks = 0:7;
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
    'L5 Centralized Request Load');

colormap(cmap);
clim([-0.5 7.5]);   % 0~7 정수값을 정확히 색 하나씩 대응
cb = colorbar;
cb.Ticks = 0:7;
cb.Label.String = 'Number of FRT Requests';

end