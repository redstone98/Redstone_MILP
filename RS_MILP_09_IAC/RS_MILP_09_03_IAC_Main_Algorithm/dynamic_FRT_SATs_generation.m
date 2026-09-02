function [rv_FRT_dynamic, FRT_ID_history, FRT_age_history, FRT_generation_history] = ... 
    dynamic_FRT_SATs_generation(number_of_FRT_SATs, tf_FRT, t_vector, rv_FRT)

%% ========================================================================
% Dynamic FRT Satellite Generation / Removal
%
% Concept:
%
%   - N FRT satellites are always maintained.
%   - Initial satellites are equally phase-spaced along one FRT period.
%   - Every phase_dt = tf_FRT / N:
%         1) Remove the leading satellite
%         2) Generate a new satellite at the t = 0 position
%
%   - Each satellite has a unique physical ID.
%
% Output:
%
%   rv_FRT_dynamic         : [nTime x 6 x N]
%   FRT_ID_history         : [nTime x N]
%   FRT_age_history        : [nTime x N]
%   FRT_generation_history : event table
% ========================================================================


%% ------------------------------------------------------------------------
% Basic parameters
% -------------------------------------------------------------------------

N_FRT = number_of_FRT_SATs;

T_FRT = tf_FRT;

phase_dt = T_FRT / N_FRT;

nTime = length(t_vector);


%% ------------------------------------------------------------------------
% Reference FRT orbit
% -------------------------------------------------------------------------
%
% rv_FRT should contain one complete FRT orbit.
%
% Ideally:
%
% size(rv_FRT,1) == length(t_vector)
%
% so that
%
% rv_FRT(k,:) = FRT state at t_vector(k)
%

if size(rv_FRT,1) == length(t_vector)

    t_FRT_reference = t_vector;

else

    % Fallback if a separate reference time vector was not saved
    t_FRT_reference = linspace( ...
        0, ...
        T_FRT, ...
        size(rv_FRT,1))';

end


%% ========================================================================
% Initial FRT constellation
% ========================================================================

initial_phase = ...
    (0:N_FRT-1)' * phase_dt;

rv0_FRT_SATs = zeros(N_FRT,6);

for sat_index = 1:N_FRT

    rv0_FRT_SATs(sat_index,:) = interp1( ...
        t_FRT_reference, ...
        rv_FRT(:,1:6), ...
        initial_phase(sat_index), ...
        'pchip');

end


%% ========================================================================
% Preallocate dynamic constellation history
% ========================================================================

rv_FRT_dynamic = NaN( ...
    nTime, ...
    6, ...
    N_FRT);

% Physical satellite ID occupying each orbital phase slot
FRT_ID_history = NaN( ...
    nTime, ...
    N_FRT);

% Age of each satellite since entering t = 0 position
FRT_age_history = NaN( ...
    nTime, ...
    N_FRT);

% Phase slot:
%
% slot 1 -> closest to t = 0 position
% slot 2 -> phase_dt ahead
% ...
% slot N -> (N-1)*phase_dt ahead


%% ========================================================================
% Dynamic constellation propagation
% ========================================================================

numerical_tol = ...
    1e-10 * max(1,T_FRT);

for time_index = 1:nTime

    current_time = t_vector(time_index);


    %% --------------------------------------------------------------------
    % Number of generation/removal events that have occurred
    %
    % event_count = 0:
    %     initial constellation
    %
    % event_count = 1:
    %     one satellite removed
    %     one new satellite generated
    %
    % etc.
    % ---------------------------------------------------------------------

    event_count = floor( ...
        (current_time + numerical_tol) / phase_dt);

    event_count = min( ...
        event_count, ...
        N_FRT);


    %% --------------------------------------------------------------------
    % Time elapsed since the latest insertion event
    % ---------------------------------------------------------------------

    time_since_event = ...
        current_time - event_count * phase_dt;

    if abs(time_since_event) < numerical_tol
        time_since_event = 0;
    end


    %% --------------------------------------------------------------------
    % Current phase of all N satellites
    %
    % Regardless of physical ID, active satellites occupy
    %
    %   tau
    %   tau + phase_dt
    %   ...
    %   tau + (N-1) phase_dt
    %
    % ---------------------------------------------------------------------

    satellite_age = ...
        time_since_event + ...
        (0:N_FRT-1)' * phase_dt;


    %% --------------------------------------------------------------------
    % Determine physical satellite IDs
    % ---------------------------------------------------------------------
    %
    % At t = 0:
    %
    % phase slot:
    %     0      dt      2dt ... (N-1)dt
    %
    % ID:
    %     1       2       3  ... N
    %
    %
    % After first event:
    %
    %     N+1     1       2  ... N-1
    %
    %
    % After second event:
    %
    %     N+2    N+1      1  ... N-2
    %
    % ---------------------------------------------------------------------

    current_IDs = zeros(N_FRT,1);

    m = event_count;

    if m == 0

        current_IDs = (1:N_FRT)';

    else

        %% Newly inserted satellites

        for slot_index = 1:m

            current_IDs(slot_index) = ...
                N_FRT + m - slot_index + 1;

        end


        %% Surviving initial satellites

        if m < N_FRT

            current_IDs(m+1:end) = ...
                (1:N_FRT-m)';

        end

    end


    %% --------------------------------------------------------------------
    % Obtain state from periodic FRT reference trajectory
    % ---------------------------------------------------------------------

    for slot_index = 1:N_FRT

        phase_time = ...
            satellite_age(slot_index);

        rv_current = interp1( ...
            t_FRT_reference, ...
            rv_FRT(:,1:6), ...
            phase_time, ...
            'pchip');

        rv_FRT_dynamic( ...
            time_index,:,slot_index) = ...
            rv_current;

    end


    %% Store ID and phase information

    FRT_ID_history(time_index,:) = ...
        current_IDs';

    FRT_age_history(time_index,:) = ...
        satellite_age';

end

%% ========================================================================
% Generation / Removal Event History
% ========================================================================

number_of_events = N_FRT;

EventNumber = (1:number_of_events)';

EventTime = ...
    EventNumber * phase_dt;

RemovedSatelliteID = zeros( ...
    number_of_events,1);

CreatedSatelliteID = zeros( ...
    number_of_events,1);


for event_index = 1:number_of_events

    %% -------------------------------------------------------------
    % Removed satellite
    %
    % Removal sequence:
    %
    % 10, 9, 8, ..., 1
    %
    % --------------------------------------------------------------

    RemovedSatelliteID(event_index) = ...
        N_FRT - event_index + 1;


    %% -------------------------------------------------------------
    % Newly created satellite
    %
    % New IDs:
    %
    % 11, 12, 13, ...
    %
    % --------------------------------------------------------------

    CreatedSatelliteID(event_index) = ...
        N_FRT + event_index;

end


FRT_generation_history = table( ...
    EventNumber, ...
    EventTime, ...
    RemovedSatelliteID, ...s
    CreatedSatelliteID);


disp(FRT_generation_history);


end