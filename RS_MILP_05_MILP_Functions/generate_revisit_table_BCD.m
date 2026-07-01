function [revisit_time_vector_info, contact_tables, revisit_vectors, satellite_cadence_info] = generate_revisit_table_BCD(access_interval_table, row_index_BCD, A_matrix, tau)


% Extract Ground Point Tables
  T = access_interval_table(row_index_BCD,:);

    % Source 열에서 고유 값 추출
    sources = unique(T.Source);
    targets = unique(T.Target);
    for i = 1:length(targets)
        src = targets(i);
        subTable = T(strcmp(T.Target, src), :);
        subTable_sorted = sortrows(subTable,4,"ascend");   
        satellite_cadence_info.(sprintf('%s_Table', src)) = subTable_sorted;
    end

    % Revisit Time Matrix (Min / Max / Mean) 생성
    Revisit_Matrix = zeros(length(sources),3);
    revisit_time_vector_combined = [];

    t0 = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');
    t_start = t0;
    t_end = t0 + seconds(max(A_matrix(:,3))+1);

    % 각각의 Source 별로 테이블을 분리하여 저장
    for i = 1:length(sources)
        src = sources{i};
        gs_index = str2double(regexp(src, '\d+', 'match', 'once'));
        % Source 값이 같은 행들만 추출
        subTable = T(strcmp(T.Source, src), :);
        subTable_sorted = sortrows(subTable,4,"ascend");

        % assignin('base', varName, subTable_sorted);
        contact_tables.(sprintf('%s_Table', src)) = subTable_sorted;

        % 기존 row에 start 및 endtime 추가
        row0 = subTable_sorted(1,:);
        
        % -------------------------
        % 윗줄 (t_start)
        % -------------------------
        row_top = row0;
        row_top.IntervalNumber = NaN;      % 필요 없으면 NaN
        row_top.StartTime = t_start;
        row_top.EndTime   = t_start;
        row_top.Duration  = 0;
        row_top.StartOrbit = 0;
        row_top.EndOrbit   = 0;
        
        % -------------------------
        % 아랫줄 (t_end)
        % -------------------------
        row_bottom = row0;
        row_bottom.IntervalNumber = NaN;
        row_bottom.StartTime = t_end;
        row_bottom.EndTime   = t_end;
        row_bottom.Duration  = 0;
        row_bottom.StartOrbit = 0;
        row_bottom.EndOrbit   = 0;
        
        % -------------------------
        % 위 + 기존 + 아래 결합
        % -------------------------
        subTable_sorted = [row_top; subTable_sorted; row_bottom];

        % Initialize Revisit Time Vector
        revisit_time_vector = zeros(height(subTable_sorted)-1,1);

        for revisit_time_index = 1:length(revisit_time_vector)
            revisit_time_vector(revisit_time_index) = seconds(table2array(subTable_sorted(revisit_time_index+1,4))-table2array(subTable_sorted(revisit_time_index,5)));

            if revisit_time_vector(revisit_time_index) < 0
               revisit_time_vector(revisit_time_index) = 0;
            end
        revisit_vectors.(sprintf('%s_revisit_time_vec', src)) = revisit_time_vector;
        end

        % 생성된 Revisit Time Vector의 Min / Max / Mean 값을 Revisit Time Matrix에 저장
        Revisit_Matrix(gs_index,1) = min(revisit_time_vector);
        Revisit_Matrix(gs_index,2) = max(revisit_time_vector);

        % 재방문 주기의 평균을 구할때는 재방문 주기가 0인 데이터 포인트를 모두 제외하였음
        revisit_time_vector = revisit_time_vector(revisit_time_vector~=0);

        Revisit_Matrix(gs_index,3) = mean(revisit_time_vector);
       revisit_time_vector_info.(['source' num2str(gs_index)]) = revisit_time_vector;
    end


    figure;
    % hold on;
    
    % for x = 1:length(Revisit_Matrix(:,1))
    %     revisit_vector = revisit_time_vector_info.(['source' num2str(x)]);
    %     n = length(revisit_vector);
    %     scatter(x*ones(n,1), revisit_vector/3600,'*','m')
    % end
   
    % Revisit Time Matrix의 값을 그래프로 출력
    x = 1:length(Revisit_Matrix(:,1));
    min_value = (Revisit_Matrix(:,1))/3600;
    max_value = (Revisit_Matrix(:,2))/3600;
    mean_value = (Revisit_Matrix(:,3))/3600;
    errorbar(x, mean_value, mean_value-min_value, max_value-mean_value,'*','LineStyle','none','color','b','MarkerEdgeColor','r')
    legend('Min, Max, Mean Revisit Time Data','Location','southoutside')    
    hold off
    title("Revisit Time Result (BCD, tau = " + tau + " sec)", 'FontSize',12,'FontWeight','bold');
    xlabel('Ground Observation Point Index','FontSize',11,'FontWeight','bold')
    ylabel('Revisit Time (Hours)','FontSize',11,'FontWeight','bold')
    xlim([-1, length(Revisit_Matrix(:,1))]+1)
        ylim([0,16])


end