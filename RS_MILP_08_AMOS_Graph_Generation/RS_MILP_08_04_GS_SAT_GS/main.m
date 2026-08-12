%% RS_MILP_08_04_GS_SAT_GS
clear; clc;
addpath ~/Desktop/Redstone_MILP/RS_MILP_08_AMOS_Graph_Generation/RS_MILP_08_04_GS_SAT_GS/
load('GS_to_SAT_Observation_Table.mat')
load('SAT_to_GS_Downlink_Table.mat')
contact_tables_BCD;
satellite_cadence_info_L2;

[imaging_to_downlink_struct, ...
 all_latency_table, ...
 ground_point_summary, ...
 overall_summary] = ...
    analyze_imaging_downlink_latency( ...
        satellite_cadence_info_L2, ...
        contact_tables_BCD);


function [imaging_to_downlink_struct, ...
          all_latency_table, ...
          ground_point_summary, ...
          overall_summary] = ...
          analyze_imaging_downlink_latency( ...
          satellite_cadence_info_L2, contact_tables_BCD)

% ANALYZE_IMAGING_DOWNLINK_LATENCY
%
% 각 위성의 촬영 시퀀스와 다운링크 컨택 시퀀스를 연결하여,
% 촬영 완료부터 가장 가까운 미래의 다운링크 가능 시점까지의
% 지연시간을 계산합니다.
%
% 입력
% -------------------------------------------------------------------------
% satellite_cadence_info_L2
%   SAT_1_Table ~ SAT_48_Table을 포함하는 struct
%   Source    : 촬영 지점, 예: "Ground_Point_1"
%   Target    : 위성, 예: "SAT_1"
%   StartTime : 촬영 시작시각
%   EndTime   : 촬영 완료시각
%
% contact_tables_BCD
%   SAT_1_Table ~ SAT_48_Table을 포함하는 struct
%   Source    : 다운링크 지상국
%   Target    : 위성
%   StartTime : 지상국 컨택 시작시각
%   EndTime   : 지상국 컨택 종료시각
%
% 출력
% -------------------------------------------------------------------------
% imaging_to_downlink_struct
%   Ground_Point_1_Table ~ Ground_Point_54_Table
%
% all_latency_table
%   모든 촬영 이벤트의 촬영-다운링크 연결 결과
%
% ground_point_summary
%   54개 촬영 지점별 이벤트 수, 다운링크 성공률 및 지연시간 통계
%
% overall_summary
%   전체 촬영 이벤트 기준 다운링크 가능 비율 및 지연시간 통계
%
% 주요 가정
% -------------------------------------------------------------------------
% 1. 촬영 영상은 동일 위성의 다운링크 컨택만 사용함.
% 2. 촬영 완료 이후 가장 먼저 이용 가능한 컨택을 선택함.
% 3. 촬영 완료시각이 컨택 시간 안에 있으면 즉시 다운링크 가능함.
% 4. 한 컨택에서 여러 영상을 모두 내려받을 수 있다고 가정함.
% 5. 데이터 용량, 전송률 및 컨택 용량 제한은 고려하지 않음.

    numSatellites   = 48;
    numGroundPoints = 54;

    % 모든 결과를 임시로 저장
    resultRows = cell(0, 10);

    % 실제로 사용된 다운링크 지상국 확인용
    downlinkStationList = strings(0,1);

    %% 위성별 촬영 이벤트와 다운링크 컨택 연결
    for satIdx = 1:numSatellites

        satField = sprintf('SAT_%d_Table', satIdx);
        satName  = sprintf("SAT_%d", satIdx);

        % 촬영 테이블이 존재하지 않는 경우
        if ~isfield(satellite_cadence_info_L2, satField)
            warning('%s가 satellite_cadence_info_L2에 없습니다.', satField);
            continue;
        end

        imagingTable = satellite_cadence_info_L2.(satField);

        if isempty(imagingTable)
            continue;
        end

        % 촬영 이벤트를 시간순으로 정렬
        imagingTable = sortrows(imagingTable, 'StartTime');

        % 다운링크 테이블 확인
        if isfield(contact_tables_BCD, satField)
            downlinkTable = contact_tables_BCD.(satField);
        else
            downlinkTable = table();
        end

        if ~isempty(downlinkTable)
            downlinkTable = sortrows(downlinkTable, 'StartTime');

            downlinkStationList = [
                downlinkStationList
                string(downlinkTable.Source)
            ];
        end

        %% 해당 위성의 각 촬영 이벤트 처리
        for imageIdx = 1:height(imagingTable)

            observationPoint = string(imagingTable.Source(imageIdx));

            imageStartTime = imagingTable.StartTime(imageIdx);
            imageEndTime   = imagingTable.EndTime(imageIdx);

            downlinkStation   = missing;
            downlinkStartTime = makeNaTLike(imageEndTime);
            downlinkEndTime   = makeNaTLike(imageEndTime);
            downlinkTime      = makeNaTLike(imageEndTime);

            latency          = seconds(NaN);
            latencySeconds   = NaN;
            isDownlinked     = false;

            if ~isempty(downlinkTable)

                % 촬영 완료 이후에도 이용 가능한 컨택 검색
                %
                % EndTime >= imageEndTime:
                % 촬영 완료시점에 이미 진행 중인 컨택도 포함
                feasibleContactIdx = find( ...
                    downlinkTable.EndTime >= imageEndTime);

                if ~isempty(feasibleContactIdx)

                    % 실제 다운링크 가능시각:
                    % max(컨택 시작시각, 촬영 완료시각)
                    availableTimes = ...
                        downlinkTable.StartTime(feasibleContactIdx);

                    activeContactIdx = availableTimes < imageEndTime;

                    % 촬영 완료 시점에 컨택이 진행 중이면
                    % 촬영 완료 즉시 다운링크 가능
                    availableTimes(activeContactIdx) = imageEndTime;

                    % 가장 빠른 이용 가능 컨택 선택
                    [downlinkTime, localIdx] = min(availableTimes);

                    selectedContactIdx = ...
                        feasibleContactIdx(localIdx);

                    downlinkStation = string( ...
                        downlinkTable.Source(selectedContactIdx));

                    downlinkStartTime = ...
                        downlinkTable.StartTime(selectedContactIdx);

                    downlinkEndTime = ...
                        downlinkTable.EndTime(selectedContactIdx);

                    latency        = downlinkTime - imageEndTime;
                    latencySeconds = seconds(latency);
                    isDownlinked   = true;
                end
            end

            resultRows(end+1,:) = { ...
                observationPoint, ...
                satName, ...
                imageStartTime, ...
                imageEndTime, ...
                downlinkStation, ...
                downlinkStartTime, ...
                downlinkEndTime, ...
                downlinkTime, ...
                latencySeconds, ...
                isDownlinked ...
                };
        end
    end

    %% 전체 결과 테이블 생성
    variableNames = { ...
        'ObservationPoint', ...
        'Satellite', ...
        'ImageStartTime', ...
        'ImageEndTime', ...
        'DownlinkStation', ...
        'ContactStartTime', ...
        'ContactEndTime', ...
        'DownlinkAvailableTime', ...
        'LatencySeconds', ...
        'IsDownlinked'};

    if isempty(resultRows)

        all_latency_table = table( ...
            strings(0,1), ...
            strings(0,1), ...
            datetime.empty(0,1), ...
            datetime.empty(0,1), ...
            strings(0,1), ...
            datetime.empty(0,1), ...
            datetime.empty(0,1), ...
            datetime.empty(0,1), ...
            zeros(0,1), ...
            false(0,1), ...
            'VariableNames', variableNames);

    else
        all_latency_table = cell2table( ...
            resultRows, ...
            'VariableNames', variableNames);
    end

    % 분석과 그래프 작성이 편리하도록 분 단위 추가
    all_latency_table.LatencyMinutes = ...
        all_latency_table.LatencySeconds / 60;

    % 시간순 정렬
    all_latency_table = sortrows( ...
        all_latency_table, ...
        {'ImageStartTime', 'Satellite'});

    %% 촬영 지점별 struct 생성
    imaging_to_downlink_struct = struct();

    for groundIdx = 1:numGroundPoints

        groundName = sprintf("Ground_Point_%d", groundIdx);
        fieldName  = sprintf('Ground_Point_%d_Table', groundIdx);

        idx = all_latency_table.ObservationPoint == groundName;

        pointTable = all_latency_table(idx,:);

        pointTable = sortrows( ...
            pointTable, ...
            {'ImageStartTime', 'Satellite'});

        imaging_to_downlink_struct.(fieldName) = pointTable;
    end

    %% 촬영 지점별 통계 계산
    pointNumber          = (1:numGroundPoints)';
    observationPointName = strings(numGroundPoints,1);

    numImagingEvents     = zeros(numGroundPoints,1);
    numDownlinked        = zeros(numGroundPoints,1);
    numNotDownlinked     = zeros(numGroundPoints,1);
    availabilityRatio    = NaN(numGroundPoints,1);

    meanLatencyMin       = NaN(numGroundPoints,1);
    medianLatencyMin     = NaN(numGroundPoints,1);
    minimumLatencyMin    = NaN(numGroundPoints,1);
    maximumLatencyMin    = NaN(numGroundPoints,1);

    for groundIdx = 1:numGroundPoints

        groundName = sprintf("Ground_Point_%d", groundIdx);
        observationPointName(groundIdx) = groundName;

        idx = all_latency_table.ObservationPoint == groundName;

        pointResult = all_latency_table(idx,:);

        numImagingEvents(groundIdx) = height(pointResult);

        if isempty(pointResult)
            continue;
        end

        successIdx = pointResult.IsDownlinked;

        numDownlinked(groundIdx) = sum(successIdx);

        numNotDownlinked(groundIdx) = ...
            numImagingEvents(groundIdx) - numDownlinked(groundIdx);

        availabilityRatio(groundIdx) = ...
            numDownlinked(groundIdx) / numImagingEvents(groundIdx);

        successfulLatency = ...
            pointResult.LatencyMinutes(successIdx);

        if ~isempty(successfulLatency)
            meanLatencyMin(groundIdx) = mean( ...
                successfulLatency, 'omitnan');

            medianLatencyMin(groundIdx) = median( ...
                successfulLatency, 'omitnan');

            minimumLatencyMin(groundIdx) = min( ...
                successfulLatency);

            maximumLatencyMin(groundIdx) = max( ...
                successfulLatency);
        end
    end

    ground_point_summary = table( ...
        pointNumber, ...
        observationPointName, ...
        numImagingEvents, ...
        numDownlinked, ...
        numNotDownlinked, ...
        availabilityRatio, ...
        availabilityRatio * 100, ...
        meanLatencyMin, ...
        medianLatencyMin, ...
        minimumLatencyMin, ...
        maximumLatencyMin, ...
        'VariableNames', { ...
        'GroundPointNumber', ...
        'ObservationPoint', ...
        'NumImagingEvents', ...
        'NumDownlinked', ...
        'NumNotDownlinked', ...
        'AvailabilityRatio', ...
        'AvailabilityPercent', ...
        'MeanLatencyMinutes', ...
        'MedianLatencyMinutes', ...
        'MinimumLatencyMinutes', ...
        'MaximumLatencyMinutes'});

    %% 전체 통계 계산
    totalImagingEvents = height(all_latency_table);

    successfulIdx = all_latency_table.IsDownlinked;

    totalDownlinked    = sum(successfulIdx);
    totalNotDownlinked = totalImagingEvents - totalDownlinked;

    successfulLatencyMinutes = ...
        all_latency_table.LatencyMinutes(successfulIdx);

    if totalImagingEvents > 0
        overallAvailabilityRatio = ...
            totalDownlinked / totalImagingEvents;
    else
        overallAvailabilityRatio = NaN;
    end

    if isempty(successfulLatencyMinutes)
        meanLatencyMinutes   = NaN;
        medianLatencyMinutes = NaN;
        minLatencyMinutes    = NaN;
        maxLatencyMinutes    = NaN;
        stdLatencyMinutes    = NaN;
    else
        meanLatencyMinutes = mean( ...
            successfulLatencyMinutes, 'omitnan');

        medianLatencyMinutes = median( ...
            successfulLatencyMinutes, 'omitnan');

        minLatencyMinutes = min(successfulLatencyMinutes);
        maxLatencyMinutes = max(successfulLatencyMinutes);

        stdLatencyMinutes = std( ...
            successfulLatencyMinutes, 'omitnan');
    end

    actualDownlinkStations = unique( ...
        downlinkStationList( ...
        ~ismissing(downlinkStationList)));

    overall_summary = struct();

    overall_summary.NumSatellites             = numSatellites;
    overall_summary.NumObservationPoints      = numGroundPoints;
    overall_summary.NumDownlinkStations       = ...
        length(actualDownlinkStations);

    overall_summary.DownlinkStationNames      = ...
        actualDownlinkStations;

    overall_summary.NumImagingEvents          = ...
        totalImagingEvents;

    overall_summary.NumDownlinkedImages       = ...
        totalDownlinked;

    overall_summary.NumNotDownlinkedImages    = ...
        totalNotDownlinked;

    overall_summary.DownlinkAvailabilityRatio = ...
        overallAvailabilityRatio;

    overall_summary.DownlinkAvailabilityPercent = ...
        overallAvailabilityRatio * 100;

    overall_summary.MeanLatencyMinutes        = ...
        meanLatencyMinutes;

    overall_summary.MedianLatencyMinutes      = ...
        medianLatencyMinutes;

    overall_summary.MinimumLatencyMinutes     = ...
        minLatencyMinutes;

    overall_summary.MaximumLatencyMinutes     = ...
        maxLatencyMinutes;

    overall_summary.StdLatencyMinutes         = ...
        stdLatencyMinutes;

    %% 결과 출력
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf(' Imaging-to-Downlink Analysis\n');
    fprintf('====================================================\n');
    fprintf('위성 수                    : %d\n', numSatellites);
    fprintf('촬영 지점 수               : %d\n', numGroundPoints);
    fprintf('확인된 다운링크 지상국 수   : %d\n', ...
        overall_summary.NumDownlinkStations);
    fprintf('전체 촬영 이벤트 수         : %d\n', ...
        totalImagingEvents);
    fprintf('다운링크 가능 영상 수        : %d\n', ...
        totalDownlinked);
    fprintf('다운링크 불가능 영상 수       : %d\n', ...
        totalNotDownlinked);
    fprintf('다운링크 가능 비율            : %.2f %%\n', ...
        overall_summary.DownlinkAvailabilityPercent);
    fprintf('평균 촬영-다운링크 시간       : %.2f min\n', ...
        meanLatencyMinutes);
    fprintf('중앙값 촬영-다운링크 시간     : %.2f min\n', ...
        medianLatencyMinutes);
    fprintf('최대 촬영-다운링크 시간       : %.2f min\n', ...
        maxLatencyMinutes);
    fprintf('====================================================\n');

    %% 전체적인 분포 그래프
    plotImagingDownlinkResults( ...
        all_latency_table, ground_point_summary);

end

%% ------------------------------------------------------------------------
function plotImagingDownlinkResults( ...
    all_latency_table, ground_point_summary)

    successIdx = all_latency_table.IsDownlinked;

    latencyMinutes = ...
        all_latency_table.LatencyMinutes(successIdx);

    if isempty(latencyMinutes)
        warning('다운링크에 성공한 촬영 이벤트가 없어 그래프를 생성하지 않습니다.');
        return;
    end

    figure( ...
        'Name', 'Imaging-to-Downlink Latency Analysis', ...
        'Position', [100 100 1250 820]);

    tiledlayout(1,2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% 1. 전체 지연시간 히스토그램
    nexttile;

    histogram( ...
        latencyMinutes, ...
        'Normalization', 'probability');

    xlabel('Imaging-to-Downlink Latency [min]');
    ylabel('Probability');
    title('Overall Latency Distribution','FontSize',12,'FontWeight','bold');

    grid on;
    box on;

    % %% 2. 경험적 누적분포함수
    % nexttile;
    % 
    % sortedLatency = sort(latencyMinutes);
    % 
    % cumulativeProbability = ...
    %     (1:length(sortedLatency))' / length(sortedLatency);
    % 
    % stairs( ...
    %     sortedLatency, ...
    %     cumulativeProbability, ...
    %     'LineWidth', 1.5);
    % 
    % xlabel('Imaging-to-Downlink Latency [min]');
    % ylabel('Cumulative Probability');
    % title('Empirical Cumulative Distributed Function (CDF) of Latency');
    % 
    % grid on;
    % box on;
    % ylim([0 1]);

    % %% 3. 촬영 지점별 다운링크 가능 비율
    % nexttile;
    % 
    % bar( ...
    %     ground_point_summary.GroundPointNumber, ...
    %     ground_point_summary.AvailabilityPercent);
    % 
    % xlabel('Observation Ground Point Number');
    % ylabel('Downlink Availability [%]');
    % title('Downlink Availability by Observation Point');
    % 
    % grid on;
    % box on;
    % 
    % ylim([0 100]);
    % xlim([0.5 54.5]);
    % xticks(1:3:54);

    %% 4. 촬영 지점별 최소-평균-최대 지연시간
    nexttile;

    groundPointNumber = ...
        ground_point_summary.GroundPointNumber;

    minLatency = ...
        ground_point_summary.MinimumLatencyMinutes;

    meanLatency = ...
        ground_point_summary.MeanLatencyMinutes;

    maxLatency = ...
        ground_point_summary.MaximumLatencyMinutes;

    % 촬영 이벤트 또는 다운링크 성공 이벤트가 없는 지점 제외
    validIdx = ...
        ~isnan(minLatency) & ...
        ~isnan(meanLatency) & ...
        ~isnan(maxLatency);

    xValue      = groundPointNumber(validIdx);
    minLatency  = minLatency(validIdx);
    meanLatency = meanLatency(validIdx);
    maxLatency  = maxLatency(validIdx);

    % 모두 열 벡터로 통일
    xValue      = xValue(:);
    minLatency  = minLatency(:);
    meanLatency = meanLatency(:);
    maxLatency  = maxLatency(:);

    % 평균값 기준 아래쪽과 위쪽 error 길이
    lowerError = meanLatency - minLatency;
    upperError = maxLatency - meanLatency;

    hold on;

    % Min-Max 범위
    hRange = errorbar( ...
        xValue, ...
        meanLatency, ...
        lowerError, ...
        upperError, ...
        'LineStyle', 'none', ...
        'LineWidth', 1.0, ...
        'CapSize', 5);

    % 평균값
    hMean = plot( ...
        xValue, ...
        meanLatency, ...
        'r.', ...
        'MarkerSize', 13);

    xlabel('Observation Ground Point Number');
    ylabel('Imaging-to-Downlink Latency [min]');

    title( ...
        'Minimum, Mean, and Maximum Latency by Observation Point','FontSize',12,'FontWeight','bold');

    grid on;
    box on;

    xlim([0.5 54.5]);
    xticks(1:3:54);

    legend( ...
        [hRange, hMean], ...
        {'Minimum-Maximum Range', 'Mean Latency'}, ...
        'Location', 'best');

    hold off;

    %% 전체 제목
    sgtitle( ...
        'Earth Observation Imaging and Downlink Sequence Analysis', ...
        'FontWeight', 'bold','fontsize',20);

end


%% ------------------------------------------------------------------------
function outputTime = makeNaTLike(referenceTime)
% 입력 datetime과 동일한 TimeZone을 갖는 NaT 생성

    if isempty(referenceTime.TimeZone)
        outputTime = NaT;
    else
        outputTime = NaT('TimeZone', referenceTime.TimeZone);
    end

end