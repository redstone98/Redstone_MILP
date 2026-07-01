%% RS_MILP_07_GS_to_SAT
function [sat_revisit_time_vector_info, satellite_contact_tables, sat_revisit_vectors, gs_cadence_info] = ...
    generate_revisit_table_ADMM(access_interval_table, row_index_ADMM, A_matrix, tau, rho)

%% Extract selected contact table
T = access_interval_table(row_index_ADMM,:);

sources = unique(T.Source);   % Ground stations
targets = unique(T.Target);   % Satellites

%% Ground station cadence info 저장
gs_cadence_info = struct();

for gs_index = 1:length(sources)
    src = sources{gs_index};
    subTable = T(strcmp(T.Source, src), :);
    subTable_sorted = sortrows(subTable, 4, "ascend");
    gs_cadence_info.(sprintf('%s_Table', src)) = subTable_sorted;
end

%% Satellite revisit matrix 생성
Satellite_Revisit_Matrix = zeros(length(targets), 3);

t0 = datetime(2030, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
t_start = t0 + seconds(min(A_matrix(:,3)));
t_end   = t0 + seconds(max(A_matrix(:,3)) + 1);

sat_revisit_time_vector_info = struct();
satellite_contact_tables = struct();
sat_revisit_vectors = struct();

%% Satellite 별로 contact sequence 계산
for sat_index = 1:length(targets)

    sat = targets{sat_index};
    sat_num = str2double(regexp(sat, '\d+', 'match', 'once'));

    subTable = T(strcmp(T.Target, sat), :);
    subTable_sorted = sortrows(subTable, 4, "ascend");

    satellite_contact_tables.(sprintf('%s_Table', sat)) = subTable_sorted;

    if isempty(subTable_sorted)
        Satellite_Revisit_Matrix(sat_index,:) = [NaN NaN NaN];
        continue;
    end

    %% Start/end boundary 추가
    row0 = subTable_sorted(1,:);

    row_top = row0;
    row_top.IntervalNumber = NaN;
    row_top.StartTime = t_start;
    row_top.EndTime   = t_start;
    row_top.Duration  = 0;
    row_top.StartOrbit = 0;
    row_top.EndOrbit   = 0;

    row_bottom = row0;
    row_bottom.IntervalNumber = NaN;
    row_bottom.StartTime = t_end;
    row_bottom.EndTime   = t_end;
    row_bottom.Duration  = 0;
    row_bottom.StartOrbit = 0;
    row_bottom.EndOrbit   = 0;

    subTable_sorted = [row_top; subTable_sorted; row_bottom];

    %% Satellite-to-any-GS revisit time 계산
    sat_revisit_time_vector = zeros(height(subTable_sorted)-1, 1);

    for k = 1:length(sat_revisit_time_vector)

        current_end_time = table2array(subTable_sorted(k, 5));
        next_start_time  = table2array(subTable_sorted(k+1, 4));

        sat_revisit_time_vector(k) = seconds(next_start_time - current_end_time);

        if sat_revisit_time_vector(k) < 0
            sat_revisit_time_vector(k) = 0;
        end
    end

    sat_revisit_vectors.(sprintf('%s_revisit_time_vec', sat)) = sat_revisit_time_vector;

    %% Min / Max / Mean
    Satellite_Revisit_Matrix(sat_num,1) = min(sat_revisit_time_vector);
    Satellite_Revisit_Matrix(sat_num,2) = max(sat_revisit_time_vector);

    sat_revisit_time_vector_nonzero = sat_revisit_time_vector(sat_revisit_time_vector ~= 0);

    if isempty(sat_revisit_time_vector_nonzero)
        Satellite_Revisit_Matrix(sat_num,3) = 0;
    else
        Satellite_Revisit_Matrix(sat_num,3) = mean(sat_revisit_time_vector_nonzero);
    end

    sat_revisit_time_vector_info.(['satellite' num2str(sat_num)]) = sat_revisit_time_vector_nonzero;

end

%% Plot
x = 1:length(targets);

min_value  = Satellite_Revisit_Matrix(:,1) / 3600;
max_value  = Satellite_Revisit_Matrix(:,2) / 3600;
mean_value = Satellite_Revisit_Matrix(:,3) / 3600;

figure;
errorbar(x, mean_value, mean_value-min_value, max_value-mean_value, ...
    '*', 'LineStyle', 'none', 'color', 'b', 'MarkerEdgeColor', 'r')

legend('Min, Max, Mean Satellite Revisit Time Data', 'Location', 'southoutside')

title(sprintf('Revisit Time Result (ADMM, tau = %d sec, \\rho = %.0e)', tau, rho), 'FontSize',12,'FontWeight','bold');

xlabel('Satellite Index', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Revisit Time to Any Ground Station (Hours)', 'FontSize', 11, 'FontWeight', 'bold')

xlim([0, length(targets)+1])

end