clear;
clc;


addpath('/Library/gurobi1300/macos_universal2/matlab')
addpath ~/Desktop/Redstone_MILP/RS_MILP_05_MILP_Functions/
addpath ~/Desktop/Redstone_MILP/RS_MILP_01_Config_MILP/
savepath
load('EOIR_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','EOIR_access_interval')

start_time = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

%% 1. Generate A matrix from contact chart
access_interval_table = EOIR_access_interval;
A_matrix = generate_A_matrix(access_interval_table, start_time);

%1.1 Unconstrained Result
[revisit_time_vector_info_uncontrained, contact_tables_uncontrained, revisit_vectors_uncontrained] = generate_revisit_table_unconstrained(access_interval_table, A_matrix, start_time);


%% 2. Input Key Parameters from A matrix


% N = number of missions
number_of_missions = 1500;
A_matrix = A_matrix(1:number_of_missions,:);
% Satellite Cadence Constraint
tau = 30;
% Number of SATs
p = 48;
% Number of GSs
q = 54;

%% 3. Generate Selection Matrices
[E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, p,q);

%% 4. Generate A, b, P for Ax=<b by given tau
[A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, p, tau, E1_Si, E2_Si_x, E2_Si_t);

%% 5. Generate E, f, g_vec for z >= Ex + f - 1
[E, f_vector, gvec] = generate_E_f_G(A_matrix, q, E1_Gj, E2_Gj_x, E2_Gj_t);

%% 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
[x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, p, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
x_BCD(abs(x_BCD) < 1e-1) = 0;
row_index_BCD = x_BCD .* A_matrix(:,4);
row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
[revisit_time_vector_info_BCD, contact_tables_BCD, revisit_vectors_BCD] = generate_revisit_table_BCD(access_interval_table, row_index_BCD, A_matrix, tau);

%% 7. Alternative Direction Multiplier Method (Use Gurobi)
maxIters   = 200;
rho        =  4 * 1e7;       % penalty parameter
[x_ADMM,z_ADMM, x_history_ADMM, z_history_ADMM] = solve_ADMM(maxIters, rho, A_matrix, p, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
row_index_ADMM = x_ADMM .* A_matrix(:,4);
row_index_ADMM = row_index_ADMM(row_index_ADMM ~= 0);
[revisit_time_vector_info_ADMM, contact_tables_ADMM, revisit_vectors_ADMM] = generate_revisit_table_ADMM(access_interval_table, row_index_ADMM, A_matrix, tau, rho);