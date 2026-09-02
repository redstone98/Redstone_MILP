function [NRHO_load_history, L4_load_history, L5_load_history] = ...
    top_1_servicing_load_history(t_vector, number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs,...
                                 best_NRHO_ID_history, best_L4_ID_history, best_L5_ID_history)

%% ========================================================================
% TOP-1 SERVICING LOAD ANALYSIS
%
% Question:
% "If every FRT independently selects its minimum-GDOP combination,
%  how many clients request the same servicing satellite?"
% ========================================================================

nTime = length(t_vector);

N_NRHO = number_of_NRHO_SATs;
N_L4   = number_of_L4_SATs;
N_L5   = number_of_L5_SATs;


%% ------------------------------------------------------------------------
% Load matrices
%
% row    = timestep
% column = servicing satellite ID
% value  = number of FRT clients requesting that satellite
% -------------------------------------------------------------------------

NRHO_load_history = zeros(nTime,N_NRHO);
L4_load_history   = zeros(nTime,N_L4);
L5_load_history   = zeros(nTime,N_L5);


for time_index = 1:nTime

    %% NRHO
    selected_NRHO = ...
        best_NRHO_ID_history(time_index,:);

    for sat_id = 1:N_NRHO

        NRHO_load_history(time_index,sat_id) = ...
            sum(selected_NRHO == sat_id);

    end


    %% L4
    selected_L4 = ...
        best_L4_ID_history(time_index,:);

    for sat_id = 1:N_L4

        L4_load_history(time_index,sat_id) = ...
            sum(selected_L4 == sat_id);

    end


    %% L5
    selected_L5 = ...
        best_L5_ID_history(time_index,:);

    for sat_id = 1:N_L5

        L5_load_history(time_index,sat_id) = ...
            sum(selected_L5 == sat_id);

    end

end


end