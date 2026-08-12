clear;
clc;


% addpath('/Library/gurobi1300/macos_universal2/matlab')
addpath ~/Desktop/Redstone_MILP/RS_MILP_06_AMOS/
addpath ~/Desktop/Redstone_MILP/RS_MILP_01_Config_MILP/
savepath
load('Scenario_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','SAT_GS_access_interval')
% load('Constellation_6G_128_SATs_16_Orbit_Planes_4_Payloads_40_inc_1hours_access_interval')
% load('EOIR_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','EOIR_access_interval')
% load('Constellation_6G_128_SATs_16_Orbit_Planes_4_Payloads_38_inc_1hours_access_interval')
start_time_original = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

%% 1. Generate A matrix from contact chart
% access_interval_table = SAT_GS_access_interval_modified;
access_interval_table = SAT_GS_access_interval;
% access_interval_table = EOIR_access_interval;
A_matrix = generate_A_matrix(access_interval_table, start_time_original);


%% 2. Input Key Parameters from A matrix
% N = number of missions

start_index = 1;
number_of_missions = 1500;
% number_of_missions = length(A_matrix(:,1));

A_matrix = A_matrix(start_index:number_of_missions,:);
access_interval_table_sorted = access_interval_table(A_matrix(:,4),:);

start_time = start_time_original + seconds(A_matrix(1,3));
end_time = start_time_original + seconds(max(A_matrix(:,3)));

t_start = A_matrix(1,3);
t_end = A_matrix(end,3);

% Satellite Cadence Constraint
tau = 4000;
% Number of SATs
number_of_SAT = 48;
% Number of GSs
number_of_GS = 54;


%2.1 Unconstrained Result
 fprintf('-----------<result>------------ \n')
 fprintf('tau = unconstrained \n');
 fprintf('total contact = %d \n' , number_of_missions);
[revisit_time_vector_info_uncontrained, contact_tables_uncontrained, revisit_vectors_uncontrained, satellite_cadence_info_unconstrained] = generate_revisit_table_unconstrained(access_interval_table_sorted, A_matrix, start_time);
 fprintf('-----------<end>-------------- \n')


%% 3. Generate Selection Matrices
[E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);

%% 4. Generate A, b, P for Ax=<b by given tau

tau_vector = tau * ones(number_of_SAT,1);
[A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, number_of_SAT, tau_vector, E1_Si, E2_Si_x, E2_Si_t);

%% 4.1 L1 Solver
% x_L1 = solve_L1(A,b,number_of_missions);
% x_L1(abs(x_L1) < 1e-1) = 0;
% row_index_L1 = x_L1 .* A_matrix(:,4);
% row_index_L1 = row_index_L1(row_index_L1 ~= 0);
% [revisit_time_vector_L1, contact_tables_L1, revisit_vectors_L1, satellite_cadence_info_L1] = generate_revisit_table_L1(access_interval_table, row_index_L1, A_matrix, tau, start_time_original);


%% 5. Generate E, f, g_vec for z >= Ex + f - 1
k_vector = ones(number_of_GS,1);
[E, f_vector, gvec_original] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end, k_vector);


% 1, 1.25, 1.5, 10
k_vector(15) = 1;

[E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_GS, E1_Gj, E2_Gj_x, E2_Gj_t, t_start, t_end, k_vector);

%% 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
[x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
z_BCD(abs(z_BCD) < 1e-1) = 0;
x_BCD(abs(x_BCD) < 1e-1) = 0;
z_BCD(abs(z_BCD) > 0.9) = 1;
row_index_BCD = x_BCD .* A_matrix(:,4);
row_index_BCD = round(row_index_BCD(row_index_BCD ~= 0));
%% 
 fprintf('-----------<Total result>------------ \n')
 fprintf('tau = %d \n',tau);
 fprintf('k_value = %4.1f \n', k_vector(15));
 fprintf('total contact = %d \n' , sum(x_BCD));
 fprintf('cost function (sec^2) = %d \n', sum(gvec_original.^2 .* z_BCD));
 fprintf('max_revisit (min) = %4.4f \n', max(gvec_original .*z_BCD)/60);
 activated_revisit_time = gvec_original .*z_BCD;
 nonzero_revisit_time = activated_revisit_time(activated_revisit_time~=0);
 fprintf('mean_revisit (min) = %4.4f \n', mean(nonzero_revisit_time/60));
 fprintf('-----------<end>-------------- \n')

%% 


[revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, satellite_cadence_info_BCD] = generate_revisit_table_BCD(access_interval_table, row_index_BCD, A_matrix, tau);



%% 7. Alternative Direction Multiplier Method (Use Gurobi)
% maxIters   = 200;
% rho        =  1 * 1e6;       % penalty parameter
% [x_ADMM,z_ADMM, x_history_ADMM, z_history_ADMM] = solve_ADMM(maxIters, rho, A_matrix, number_of_SAT, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info);
% row_index_ADMM = x_ADMM .* A_matrix(:,4);
% row_index_ADMM = row_index_ADMM(row_index_ADMM ~= 0);
% [revisit_time_vector_info_ADMM, contact_tables_ADMM, revisit_vectors_ADMM, satellite_cadence_info_ADMM] = generate_revisit_table_ADMM(access_interval_table, row_index_ADMM, A_matrix, tau, rho);



%% 8. Plot the Revisit Block Graph for single GS


%% 8.1 Max |x|
% figure; hold on; grid on;
% 
% T = contact_tables_L1.Ground_Point_15_Table;
% 
% % StartTime, EndTime이 datetime이라고 가정
% t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
% t_mid = [start_time; t_mid; end_time];
% 
% % 시간순 정렬
% t_mid = sort(t_mid);
% 
% % 다음 contact까지의 시간 차이 [sec]
% delta_t = seconds(diff(t_mid));
% 
% % 마지막 점은 다음 점이 없으므로 제외
% t_plot = t_mid(1:end-1);
% 
% for k = 1:length(delta_t)
% 
%     % x축 폭도 delta_t 만큼
%     x0 = t_plot(k);
%     x1 = t_plot(k) + seconds(delta_t(k));
% 
%     % y축 높이도 delta_t 만큼
%     y0 = 0;
%     y1 = delta_t(k);
% 
%     patch([x0 x1 x1 x0], ...
%           [y0 y0 y1 y1]/60, ...
%           0.8*ones(1,3), ...
%           'FaceColor', 'b', ...
%           'FaceAlpha', 0.1, ...
%           'EdgeColor', 'r', ...
%           'LineWidth',1.5);
% end
% 
% xlabel('Time','FontSize',11,'FontWeight','bold');
% ylabel('\Delta t to next contact [min]','FontSize',11,'FontWeight','bold');
% % ylim([0,700])
% title('[Max |x|] Revisit time for Ground Point 10','FontSize',12,'FontWeight','bold');

%% 8.3 BCD

figure; hold on; grid on;
T = contact_tables_BCD.Ground_Point_15_Table;

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
xlabel('Time','FontSize',11,'FontWeight','bold');
ylabel('\Delta t to next contact [min]','FontSize',11,'FontWeight','bold');
% ylim([0,700])
title('[BCD] Revisit time for Ground Point 15','FontSize',12,'FontWeight','bold');

 delta_t = delta_t(delta_t ~= 0);
 fprintf('-----------<GS15 result>------------ \n')
 fprintf('tau = %d \n',tau);
 fprintf('k value = %4.4f \n', k_vector(15));
 fprintf('GS15 total contact = %d \n' , length(delta_t)-1);
 fprintf('GS15 cost function (sec^2) = %d \n', sum(delta_t.^2));
 fprintf('GS15 max_revisit (min) = %4.4f \n', max(delta_t)/60);
 fprintf('GS15 mean_revisit (min) = %4.4f \n', mean(delta_t)/60);
 fprintf('-----------<end>-------------- \n')

%% 8.4 Unconstrained


figure; hold on; grid on;
T = contact_tables_uncontrained.Ground_Point_15_Table;

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

xlabel('Time','FontSize',11,'FontWeight','bold');
ylabel('\Delta t to next contact [min]','FontSize',11,'FontWeight','bold');
% ylim([0,700])
title('[Unconstrained]Revisit time for Ground Point 10','FontSize',12,'FontWeight','bold');

 delta_t = delta_t(delta_t ~= 0);
 fprintf('-----------<GS15 result>------------ \n')
 fprintf('tau = unconstrained \n');
 fprintf('GS15 total contact = %d \n' , length(delta_t)-1);
 fprintf('GS15 cost function (sec^2) = %d \n', sum(delta_t.^2));
 fprintf('GS15 max_revisit (min) = %4.4f \n', max(delta_t)/60);
 fprintf('GS15 mean_revisit (min) = %4.4f \n', mean(delta_t)/60);
 fprintf('-----------<end>-------------- \n')
