%% RS_MILP_07_GS_to_SAT

clear;
clc;


% addpath('/Library/gurobi1300/macos_universal2/matlab')
addpath ~/Desktop/Redstone_MILP/RS_MILP_07_GS_to_SAT/
savepath
load('Scenario_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','SAT_GS_access_interval')

start_time_original = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

%% 1. Generate A matrix from contact chart
access_interval_table = SAT_GS_access_interval;
A_matrix = generate_A_matrix(access_interval_table, start_time_original);


%% 2. Input Key Parameters from A matrix
% N = number of missions

start_index = 1;
number_of_missions = 1500;

A_matrix = A_matrix(start_index:number_of_missions,:);
access_interval_table_sorted = access_interval_table(A_matrix(:,4),:);

start_time = start_time_original + seconds(A_matrix(1,3));
end_time = start_time_original + seconds(max(A_matrix(:,3)));

t_start = A_matrix(1,3);
t_end = A_matrix(end,3);

% Satellite Cadence Constraint
tau = 300;
% Number of SATs
number_of_SAT = 48;
% Number of GSs
number_of_GS = 31;


%2.1 Unconstrained Result
[revisit_time_vector_info_uncontrained, contact_tables_uncontrained, revisit_vectors_uncontrained, satellite_cadence_info_unconstrained] = generate_revisit_table_unconstrained(access_interval_table_sorted, A_matrix, start_time);



%% 3. Generate Selection Matrices
[E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);

%% 4. Generate A, b, P for Ax=<b by given tau

tau_vector = zeros(number_of_GS,1);
tau_vector(1:31) = 3000;


[A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);

%% 5. Generate E, f, g_vec for z >= Ex + f - 1
k_vector = ones(number_of_SAT,1) * 0.1;
k_vector(10) = 0.1;
k_vector(30) = 0.1;


[E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);

%% 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
[x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
x_BCD(abs(x_BCD) < 1e-1) = 0;
row_index_BCD = x_BCD .* A_matrix(:,4);
row_index_BCD = row_index_BCD(row_index_BCD ~= 0);


[revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD(access_interval_table, row_index_BCD, A_matrix, tau);



%% 7. Alternative Direction Multiplier Method (Use Gurobi)
maxIters   = 200;
rho        =  1 * 1e6;       % penalty parameter
[x_ADMM,z_ADMM,x_hist,z_hist] = solve_ADMM(maxIters, rho, A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec, E1_Gj, S_j_info, U_j_info);
row_index_ADMM = x_ADMM .* A_matrix(:,4);
row_index_ADMM = row_index_ADMM(row_index_ADMM ~= 0);
[revisit_time_vector_info_ADMM, contact_tables_ADMM, revisit_vectors_ADMM, GS_cadence_info_ADMM] = generate_revisit_table_ADMM(access_interval_table, row_index_ADMM, A_matrix, tau, rho);



%% 8. Plot the Revisit Block Graph for single GS
figure; hold on; grid on;

T = contact_tables_ADMM.SAT_23_Table;

% StartTime, EndTime이 datetime이라고 가정
t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
t_mid = [start_time; t_mid; end_time];

% 시간순 정렬
t_mid = sort(t_mid);

% 다음 contact까지의 시간 차이 [sec]
delta_t = seconds(diff(t_mid));

% 마지막 점은 다음 점이 없으므로 제외
t_plot = t_mid(1:end-1);

for k = 1:length(delta_t)

    % x축 폭도 delta_t 만큼
    x0 = t_plot(k);
    x1 = t_plot(k) + seconds(delta_t(k));

    % y축 높이도 delta_t 만큼
    y0 = 0;
    y1 = delta_t(k);

    patch([x0 x1 x1 x0], ...
          [y0 y0 y1 y1]/60, ...
          0.8*ones(1,3), ...
          'FaceColor', 'c', ...
          'FaceAlpha', 0.1, ...
          'EdgeColor', 'g', ...
          'LineWidth',1.5);
end

xlabel('Time');
ylabel('\Delta t to next contact [min]');
title('t-\Delta t Square Plot');


T = contact_tables_BCD.SAT_23_Table;

% StartTime, EndTime이 datetime이라고 가정
t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
t_mid = [start_time; t_mid; end_time];

% 시간순 정렬
t_mid = sort(t_mid);

% 다음 contact까지의 시간 차이 [sec]
delta_t = seconds(diff(t_mid));

% 마지막 점은 다음 점이 없으므로 제외
t_plot = t_mid(1:end-1);

for k = 1:length(delta_t)

    % x축 폭도 delta_t 만큼
    x0 = t_plot(k);
    x1 = t_plot(k) + seconds(delta_t(k));

    % y축 높이도 delta_t 만큼
    y0 = 0;
    y1 = delta_t(k);

    patch([x0 x1 x1 x0], ...
          [y0 y0 y1 y1]/60, ...
          0.8*ones(1,3), ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.1, ...
          'EdgeColor', 'r', ...
          'LineWidth',1.5);
end

T = contact_tables_uncontrained.SAT_23_Table;

% StartTime, EndTime이 datetime이라고 가정
t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
t_mid = [start_time; t_mid; end_time];

% 시간순 정렬
t_mid = sort(t_mid);

% 다음 contact까지의 시간 차이 [sec]
delta_t = seconds(diff(t_mid));

% 마지막 점은 다음 점이 없으므로 제외
t_plot = t_mid(1:end-1);

for k = 1:length(delta_t)

    % x축 폭도 delta_t 만큼
    x0 = t_plot(k);
    x1 = t_plot(k) + seconds(delta_t(k));

    % y축 높이도 delta_t 만큼
    y0 = 0;
    y1 = delta_t(k);

    patch([x0 x1 x1 x0], ...
          [y0 y0 y1 y1]/60, ...
          0.8*ones(1,3), ...
          'FaceColor', 'r', ...
          'FaceAlpha', 0.1, ...
          'EdgeColor', 'b', ...
          'LineWidth',1);
end
