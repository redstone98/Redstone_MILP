%% RS_MILP_08_03_GS_to_SAT

clear;
clc;


% addpath('/Library/gurobi1300/macos_universal2/matlab')
addpath ~/Desktop/Redstone_MILP/RS_MILP_08_AMOS_Graph_Generation/RS_MILP_08_03_GS_to_SAT/
savepath
load('Scenario_48_SATs_9_GS_12_Orbit_Planes_98_inc_1days_access_interval.mat','SAT_GS_access_interval')

start_time_original = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

%% 1. Generate A matrix from contact chart
access_interval_table = SAT_GS_access_interval;
A_matrix = generate_A_matrix(access_interval_table, start_time_original);


%% 2. Input Key Parameters from A matrix
% N = number of missions

start_index = 1;
number_of_missions = length(A_matrix(:,1));

A_matrix = A_matrix(start_index:number_of_missions,:);
access_interval_table_sorted = access_interval_table(A_matrix(:,4),:);

start_time = start_time_original + seconds(A_matrix(1,3));
end_time = start_time_original + seconds(max(A_matrix(:,3)));

t_start = A_matrix(1,3);
t_end = A_matrix(end,3);

% Cadence Constraint
tau = 1500;
% Numbe of SATs
number_of_SAT = 48;
% Number of GSs
number_of_GS = 9;


%2.1 Unconstrained Result
 fprintf('-----------<result>------------ \n')
 fprintf('tau = unconstrained \n');
 fprintf('total contact = %d \n' , number_of_missions);
[revisit_time_vector_info_uncontrained, contact_tables_uncontrained, revisit_vectors_uncontrained, GS_cadence_info_unconstrained] = generate_revisit_table_unconstrained(access_interval_table_sorted, A_matrix, start_time);
 fprintf('-----------<end>-------------- \n')


%% 3. Generate Selection Matrices
[E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, number_of_SAT,number_of_GS);



%% tau_selected = 1500 sec

% 4. Generate A, b, P for Ax=<b by given tau
tau_selected = tau;
tau_vector = zeros(number_of_GS,1);
tau_vector(1:number_of_GS) = tau;
% tau_vector(2) = tau_selected;

[A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);

% 5. Generate E, f, g_vec for z >= Ex + f - 1
k_vector = ones(number_of_SAT,1);
[E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);

% 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
[x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
z_BCD(abs(z_BCD) < 1e-1) = 0;
x_BCD(abs(x_BCD) < 1e-1) = 0;
z_BCD(abs(z_BCD) > 0.9) = 1;
row_index_BCD = x_BCD .* A_matrix(:,4);
row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
row_index_BCD = round(row_index_BCD);


 fprintf('-----------<Total result>------------ \n')
 fprintf('tau = %d \n',tau);
 fprintf('total contact = %d \n' , sum(x_BCD));
 fprintf('cost function (sec^2) = %d \n', sum(gvec.^2 .* z_BCD));
 fprintf('max_revisit (min) = %4.4f \n', max(gvec .*z_BCD)/60);
 activated_revisit_time = gvec .*z_BCD;
 nonzero_revisit_time = activated_revisit_time(activated_revisit_time~=0);
 fprintf('mean_revisit (min) = %4.4f \n', mean(nonzero_revisit_time/60));
 fprintf('-----------<end>-------------- \n')


[revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD_selected(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);
% [revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);

% 8. Plot the Revisit Block Graph for single GS
figure; hold on; grid on;
T = GS_cadence_info_BCD.Ground_Point_2_Table;

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
          [y0 y0 y1 y1], ...
          0.8*ones(1,3), ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.1, ...
          'EdgeColor', 'r', ...
          'LineWidth',1.5);
end
xlabel('Time','FontSize',11,'FontWeight','bold');
ylabel('\Delta t to next contact [seconds]','FontSize',11,'FontWeight','bold');
xlim([t_plot(1), t_plot(k) + seconds(delta_t(k))])
ylim([0,4100])
title("Operation Cadence for GS 2 [Svalbard, tau_2 = "+ tau_selected +" secs]","FontSize",12,"FontWeight","bold");


 delta_t = delta_t(delta_t ~= 0);
 fprintf('-----------<Svarbard Single GS result>------------ \n')
 fprintf('total contact = %d \n' , length(delta_t)-1);
 fprintf('cost function (sec^2) = %d \n', sum(delta_t.^2));
 fprintf('mean_cadence (min) = %4.4f \n', mean(delta_t)/60);
  fprintf('max_cadence (min) = %4.4f \n', max(delta_t)/60);
 fprintf('-----------<end>-------------- \n')

% %% tau_selected = 1400 sec
% 
% % 4. Generate A, b, P for Ax=<b by given tau
% tau_selected = 1400;
% tau_vector = zeros(number_of_GS,1);
% tau_vector(1:number_of_GS) = 1500;
% tau_vector(2) = tau_selected;
% 
% [A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);
% 
% % 5. Generate E, f, g_vec for z >= Ex + f - 1
% k_vector = ones(number_of_SAT,1);
% [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);
% 
% % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% [x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
% x_BCD(abs(x_BCD) < 1e-1) = 0;
% row_index_BCD = x_BCD .* A_matrix(:,4);
% row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
% row_index_BCD = round(row_index_BCD);
% 
% 
% 
% [revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD_selected(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);
% % 8. Plot the Revisit Block Graph for single GS
% figure; hold on; grid on;
% T = GS_cadence_info_BCD.Ground_Point_2_Table;
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
%           [y0 y0 y1 y1], ...
%           0.8*ones(1,3), ...
%           'FaceColor', 'b', ...
%           'FaceAlpha', 0.1, ...
%           'EdgeColor', 'r', ...
%           'LineWidth',1.5);
% end
% xlabel('Time','FontSize',11,'FontWeight','bold');
% ylabel('\Delta t to next contact [seconds]','FontSize',11,'FontWeight','bold');
% xlim([t_plot(1), t_plot(k) + seconds(delta_t(k))])
% ylim([0,3500])
% title("Operation Cadence for GS 2 [Svalbard, tau_2 = "+ tau_selected +" secs]","FontSize",12,"FontWeight","bold");
% 
% % %% tau_selected = 1300 sec
% % 
% % % 4. Generate A, b, P for Ax=<b by given tau
% % tau_selected = 1300;
% % tau_vector = zeros(number_of_GS,1);
% % tau_vector(1:number_of_GS) = 1500;
% % tau_vector(2) = tau_selected;
% % 
% % [A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);
% % 
% % % 5. Generate E, f, g_vec for z >= Ex + f - 1
% % k_vector = ones(number_of_SAT,1);
% % [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);
% % 
% % % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% % [x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
% % x_BCD(abs(x_BCD) < 1e-1) = 0;
% % row_index_BCD = x_BCD .* A_matrix(:,4);
% % row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
% % row_index_BCD = round(row_index_BCD);
% % 
% % [revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD_selected(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);% 8. Plot the Revisit Block Graph for single GS
% % figure; hold on; grid on;
% % T = GS_cadence_info_BCD.Ground_Point_2_Table;
% % 
% % % StartTime, EndTime이 datetime이라고 가정
% % t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
% % t_mid = [start_time; t_mid; end_time];
% % 
% % % 시간순 정렬
% % t_mid = sort(t_mid);
% % 
% % % 다음 contact까지의 시간 차이 [sec]
% % delta_t = seconds(diff(t_mid));
% % 
% % % 마지막 점은 다음 점이 없으므로 제외
% % t_plot = t_mid(1:end-1);
% % 
% % for k = 1:length(delta_t)
% % 
% %     % x축 폭도 delta_t 만큼
% %     x0 = t_plot(k);
% %     x1 = t_plot(k) + seconds(delta_t(k));
% % 
% %     % y축 높이도 delta_t 만큼
% %     y0 = 0;
% %     y1 = delta_t(k);
% % 
% %     patch([x0 x1 x1 x0], ...
% %           [y0 y0 y1 y1], ...
% %           0.8*ones(1,3), ...
% %           'FaceColor', 'b', ...
% %           'FaceAlpha', 0.1, ...
% %           'EdgeColor', 'r', ...
% %           'LineWidth',1.5);
% % end
% % xlabel('Time','FontSize',11,'FontWeight','bold');
% % ylabel('\Delta t to next contact [seconds]','FontSize',11,'FontWeight','bold');
% % xlim([t_plot(1), t_plot(k) + seconds(delta_t(k))])
% % ylim([0,3500])
% % title("Operation Cadence for GS 2 [Svalbard, tau_2 = "+ tau_selected +" secs]","FontSize",12,"FontWeight","bold");
% % 
% % %% tau_selected = 1200 sec
% % 
% % % 4. Generate A, b, P for Ax=<b by given tau
% % tau_selected = 1200;
% % tau_vector = zeros(number_of_GS,1);
% % tau_vector(1:number_of_GS) = 1500;
% % tau_vector(2) = tau_selected;
% % 
% % [A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);
% % 
% % % 5. Generate E, f, g_vec for z >= Ex + f - 1
% % k_vector = ones(number_of_SAT,1);
% % [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);
% % 
% % % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% % [x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
% % x_BCD(abs(x_BCD) < 1e-1) = 0;
% % row_index_BCD = x_BCD .* A_matrix(:,4);
% % row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
% % row_index_BCD = round(row_index_BCD);
% % 
% % [revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD_selected(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);
% % % 8. Plot the Revisit Block Graph for single GS
% % figure; hold on; grid on;
% % T = GS_cadence_info_BCD.Ground_Point_2_Table;
% % 
% % % StartTime, EndTime이 datetime이라고 가정
% % t_mid = T.StartTime + (T.EndTime - T.StartTime)/2;
% % t_mid = [start_time; t_mid; end_time];
% % 
% % % 시간순 정렬
% % t_mid = sort(t_mid);
% % 
% % % 다음 contact까지의 시간 차이 [sec]
% % delta_t = seconds(diff(t_mid));
% % 
% % % 마지막 점은 다음 점이 없으므로 제외
% % t_plot = t_mid(1:end-1);
% % 
% % for k = 1:length(delta_t)
% % 
% %     % x축 폭도 delta_t 만큼
% %     x0 = t_plot(k);
% %     x1 = t_plot(k) + seconds(delta_t(k));
% % 
% %     % y축 높이도 delta_t 만큼
% %     y0 = 0;
% %     y1 = delta_t(k);
% % 
% %     patch([x0 x1 x1 x0], ...
% %           [y0 y0 y1 y1], ...
% %           0.8*ones(1,3), ...
% %           'FaceColor', 'b', ...
% %           'FaceAlpha', 0.1, ...
% %           'EdgeColor', 'r', ...
% %           'LineWidth',1.5);
% % end
% % xlabel('Time','FontSize',11,'FontWeight','bold');
% % ylabel('\Delta t to next contact [seconds]','FontSize',11,'FontWeight','bold');
% % xlim([t_plot(1), t_plot(k) + seconds(delta_t(k))])
% % ylim([0,3500])
% % title("Operation Cadence for GS 2 [Svalbard, tau_2 = "+ tau_selected +" secs]","FontSize",12,"FontWeight","bold");
% 
% 
% %% tau_selected = 1100 sec
% 
% % 4. Generate A, b, P for Ax=<b by given tau
% tau_selected = 1100;
% tau_vector = zeros(number_of_GS,1);
% tau_vector(1:number_of_GS) = 1500;
% tau_vector(2) = tau_selected;
% 
% [A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, number_of_GS, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t);
% 
% % 5. Generate E, f, g_vec for z >= Ex + f - 1
% k_vector = ones(number_of_SAT,1);
% [E, f_vector, gvec] = generate_E_f_G(A_matrix, number_of_SAT, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector);
% 
% % 6. Block Coordiate Decent Optimization Solver (Use HiGHS)
% [x_BCD,z_BCD, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, number_of_GS, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info);
% x_BCD(abs(x_BCD) < 1e-1) = 0;
% row_index_BCD = x_BCD .* A_matrix(:,4);
% row_index_BCD = row_index_BCD(row_index_BCD ~= 0);
% row_index_BCD = round(row_index_BCD);
% 
% [revisit_time_vector_info, contact_tables_BCD, revisit_vectors_BCD, GS_cadence_info_BCD] = generate_revisit_table_BCD_selected(access_interval_table, row_index_BCD, A_matrix, tau, tau_selected);
% % 8. Plot the Revisit Block Graph for single GS
% figure; hold on; grid on;
% T = GS_cadence_info_BCD.Ground_Point_2_Table;
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
%           [y0 y0 y1 y1], ...
%           0.8*ones(1,3), ...
%           'FaceColor', 'b', ...
%           'FaceAlpha', 0.1, ...
%           'EdgeColor', 'r', ...
%           'LineWidth',1.5);
% end
% xlabel('Time','FontSize',11,'FontWeight','bold');
% ylabel('\Delta t to next contact [seconds]','FontSize',11,'FontWeight','bold');
% xlim([t_plot(1), t_plot(k) + seconds(delta_t(k))])
% ylim([0,4100])
% title("Operation Cadence for GS 2 [Svalbard, tau_2 = "+ tau_selected +" secs]","FontSize",12,"FontWeight","bold");