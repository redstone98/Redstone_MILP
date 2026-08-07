clear;
clc;


% addpath('/Library/gurobi1300/macos_universal2/matlab')
addpath ~/Desktop/Redstone_MILP/RS_MILP_08_AMOS_Graph_Generation/RS_MILP_08_02_SAT_to_GS/
addpath ~/Desktop/Redstone_MILP/RS_MILP_01_Config_MILP/
addpath ~/Desktop/Redstone_MILP/RS_MILP_06_AMOS/
savepath = '~/Desktop/Redstone_MILP/RS_MILP_08_AMOS_Graph_Generation/RS_MILP_08_02_SAT_to_GS/';
load('Scenario_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','SAT_GS_access_interval')
% load('EOIR_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','EOIR_access_interval')

start_time_original = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

%% 1. Generate A matrix from contact chart
access_interval_table = SAT_GS_access_interval;
% access_interval_table = EOIR_access_interval;
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


% Number of SATs
number_of_SAT = 48;
% Number of GSs
number_of_GS = 54;


%2.1 Unconstrained Result
[revisit_time_vector_info_uncontrained, contact_tables_uncontrained, revisit_vectors_uncontrained, satellite_cadence_info_unconstrained] = generate_revisit_table_unconstrained(access_interval_table_sorted, A_matrix, start_time_original);


%% Satellite Cadence Constraint = 1000
tau = 1000;
% 3. Generate Selection Matrices
[E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);

% 4. Generate A, b, P for Ax=<b by given tau
[A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, number_of_SAT, tau, E1_Si, E2_Si_x, E2_Si_t);

% 5. Generate E, f, g_vec for z >= Ex + f - 1
[E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end);

% 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
[x_L2,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
x_L2(abs(x_L2) < 1e-1) = 0;
row_index_L2 = x_L2 .* A_matrix(:,4);
row_index_L2 = row_index_L2(row_index_L2 ~= 0);
[revisit_time_vector_L2, contact_tables_L2, revisit_vectors_L2, satellite_cadence_info_L2] = generate_revisit_table_L2(access_interval_table, row_index_L2, A_matrix, tau);

% 
% filename = 'GS_to_SAT_Observation_Table';
% fullname = fullfile(savepath,filename);
% save(fullname,'satellite_cadence_info_L2')


% %% Satellite Cadence Constraint = 2000
% tau = 2000;
% % 3. Generate Selection Matrices
% [E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);
% 
% % 4. Generate A, b, P for Ax=<b by given tau
% [A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, number_of_SAT, tau, E1_Si, E2_Si_x, E2_Si_t);
% 
% % 5. Generate E, f, g_vec for z >= Ex + f - 1
% [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end);
% 
% % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% [x_L2,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
% x_L2(abs(x_L2) < 1e-1) = 0;
% row_index_L2 = x_L2 .* A_matrix(:,4);
% row_index_L2 = row_index_L2(row_index_L2 ~= 0);
% [revisit_time_vector_L2, contact_tables_L2, revisit_vectors_L2, satellite_cadence_info_L2] = generate_revisit_table_L2(access_interval_table, row_index_L2, A_matrix, tau);
% 
% %% Satellite Cadence Constraint = 3000
% tau = 3000;
% % 3. Generate Selection Matrices
% [E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);
% 
% % 4. Generate A, b, P for Ax=<b by given tau
% [A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, number_of_SAT, tau, E1_Si, E2_Si_x, E2_Si_t);
% 
% % 5. Generate E, f, g_vec for z >= Ex + f - 1
% [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end);
% 
% % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% [x_L2,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
% x_L2(abs(x_L2) < 1e-1) = 0;
% row_index_L2 = x_L2 .* A_matrix(:,4);
% row_index_L2 = row_index_L2(row_index_L2 ~= 0);
% [revisit_time_vector_L2, contact_tables_L2, revisit_vectors_L2, satellite_cadence_info_L2] = generate_revisit_table_L2(access_interval_table, row_index_L2, A_matrix, tau);
% 
% %% Satellite Cadence Constraint = 4000
% tau = 4000;
% % 3. Generate Selection Matrices
% [E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);
% 
% % 4. Generate A, b, P for Ax=<b by given tau
% [A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, number_of_SAT, tau, E1_Si, E2_Si_x, E2_Si_t);
% 
% % 5. Generate E, f, g_vec for z >= Ex + f - 1
% [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end);
% 
% % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% [x_L2,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
% x_L2(abs(x_L2) < 1e-1) = 0;
% row_index_L2 = x_L2 .* A_matrix(:,4);
% row_index_L2 = row_index_L2(row_index_L2 ~= 0);
% [revisit_time_vector_L2, contact_tables_L2, revisit_vectors_L2, satellite_cadence_info_L2] = generate_revisit_table_L2(access_interval_table, row_index_L2, A_matrix, tau);