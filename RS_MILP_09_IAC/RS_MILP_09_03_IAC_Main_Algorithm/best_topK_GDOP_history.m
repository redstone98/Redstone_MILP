function [FRT_GDOP_lookup_table, ...
 best_GDOP_history, best_combination_index_history, best_NRHO_ID_history, best_L4_ID_history, best_L5_ID_history, ...
 topK_GDOP_history, topK_combination_index_history, topK_NRHO_ID_history, topK_L4_ID_history, topK_L5_ID_history] = ...
 best_topK_GDOP_history(TOP_K, t_vector, tf_FRT, number_of_FRT_SATs, number_of_NRHO_SATs, number_of_L4_SATs, number_of_L5_SATs, ...
                        rv_NRHO_constellation, rv_L4_constellation, rv_L5_constellation, FRT_ID_history, FRT_age_history, rv_FRT_dynamic)




%% ========================================================================
% ID-AWARE OPTIMAL GDOP LOOKUP FOR ALL ACTIVE FRT SATELLITES
%
% For every timestep and every active FRT satellite:
%
%   FRT_i -> Earth + NRHO_j + L4_k + L5_l
%
% Find:
%
%   minimum GDOP
%   best NRHO ID
%   best L4 ID
%   best L5 ID
%   best combination index
%
% Also retain Top-K candidates for later task allocation.
% ========================================================================
global mu

%% ------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------

nTime = length(t_vector);

N_FRT  = number_of_FRT_SATs;
N_NRHO = number_of_NRHO_SATs;
N_L4   = number_of_L4_SATs;
N_L5   = number_of_L5_SATs;

nCombination = ...
    N_NRHO * N_L4 * N_L5;

% Numerical singularity tolerance
GDOP_rcond_tolerance = 1e-12;


%% ------------------------------------------------------------------------
% Earth position
% -------------------------------------------------------------------------
%
% IMPORTANT:
%
% If your CRTBP frame is Earth-centered:
%       r_Earth = [0;0;0]
%
% If it is the standard barycentric CR3BP rotating frame:
%       r_Earth = [-mu;0;0]
%

r_Earth = [-mu; 0; 0];


%% ========================================================================
% Combination index mapping
%
% CombinationIndex -> NRHO / L4 / L5
% ========================================================================

comb_NRHO = zeros(nCombination,1);
comb_L4   = zeros(nCombination,1);
comb_L5   = zeros(nCombination,1);

combination_index = 0;

for NRHO_index = 1:N_NRHO
    for L4_index = 1:N_L4
        for L5_index = 1:N_L5

            combination_index = combination_index + 1;

            comb_NRHO(combination_index) = NRHO_index;
            comb_L4(combination_index)   = L4_index;
            comb_L5(combination_index)   = L5_index;

        end
    end
end


%% ========================================================================
% Preallocate BEST solution history
% ========================================================================

best_GDOP_history = NaN(nTime,N_FRT);

best_combination_index_history = NaN(nTime,N_FRT);

best_NRHO_ID_history = NaN(nTime,N_FRT);
best_L4_ID_history   = NaN(nTime,N_FRT);
best_L5_ID_history   = NaN(nTime,N_FRT);


%% ========================================================================
% Preallocate Top-K history
%
% dimensions:
%
%   time x FRT slot x candidate rank
% ========================================================================

topK_GDOP_history = NaN( ...
    nTime, ...
    N_FRT, ...
    TOP_K);

topK_combination_index_history = NaN( ...
    nTime, ...
    N_FRT, ...
    TOP_K);

topK_NRHO_ID_history = NaN( ...
    nTime, ...
    N_FRT, ...
    TOP_K);

topK_L4_ID_history = NaN( ...
    nTime, ...
    N_FRT, ...
    TOP_K);

topK_L5_ID_history = NaN( ...
    nTime, ...
    N_FRT, ...
    TOP_K);


%% ========================================================================
% MAIN GDOP LOOP
% ========================================================================

I4 = eye(4);

for time_index = 1:nTime

    fprintf( ...
        'GDOP calculation: timestep %d / %d\n', ...
        time_index, ...
        nTime);


    %% --------------------------------------------------------------------
    % Servicing satellite positions at current timestep
    % ---------------------------------------------------------------------

    r_NRHO_all = squeeze( ...
        rv_NRHO_constellation(time_index,1:3,:))';

    r_L4_all = squeeze( ...
        rv_L4_constellation(time_index,1:3,:))';

    r_L5_all = squeeze( ...
        rv_L5_constellation(time_index,1:3,:))';


    %% ====================================================================
    % Loop through all currently active FRT satellites
    % ====================================================================

    for FRT_slot = 1:N_FRT

        %% ================================================================
        % Current FRT satellite
        % ================================================================
        
        current_FRT_ID = ...
            FRT_ID_history(time_index,FRT_slot);
        
        
        %% -------------------------------------------------------------
        % Current FRT position
        %
        % IMPORTANT:
        % Always force it to 1 x 3 row vector
        % -------------------------------------------------------------
        
        r_FRT = reshape( ...
            rv_FRT_dynamic(time_index,1:3,FRT_slot), ...
            1,3);
        
        
        %% ================================================================
        % Earth LOS
        % ================================================================
        
        % Earth = 3 x 1 column
        % r_FRT' = 3 x 1 column
        
        rho_Earth = ...
            r_Earth - r_FRT';
        
        d_Earth = ...
            norm(rho_Earth);
        
        if d_Earth <= eps
            continue;
        end
        
        u_Earth = ...
            rho_Earth / d_Earth;
        
        % Force 3 x 1
        u_Earth = ...
            u_Earth(:);
        
        % 4 x 1
        h_Earth = [
            u_Earth
            1
            ];
        
        % 4 x 4
        M_Earth = ...
            h_Earth * h_Earth';
        
        
        %% ================================================================
        % NRHO LOS
        % ================================================================
        %
        % r_NRHO_all : N_NRHO x 3
        % r_FRT      : 1 x 3
        %
        % MATLAB implicit expansion:
        %
        % [N_NRHO x 3] - [1 x 3]
        %
        % => [N_NRHO x 3]
        %
        
        rho_NRHO = ...
            r_NRHO_all - r_FRT;
        
        d_NRHO = ...
            vecnorm( ...
                rho_NRHO, ...
                2, ...
                2);
        
        u_NRHO = ...
            rho_NRHO ./ d_NRHO;
        
        H_NRHO = [
            u_NRHO, ...
            ones(N_NRHO,1)
            ];
        
        
        %% ================================================================
        % L4 LOS
        % ================================================================
        
        rho_L4 = ...
            r_L4_all - r_FRT;
        
        d_L4 = ...
            vecnorm( ...
                rho_L4, ...
                2, ...
                2);
        
        u_L4 = ...
            rho_L4 ./ d_L4;
        
        H_L4 = [
            u_L4, ...
            ones(N_L4,1)
            ];
        
        
        %% ================================================================
        % L5 LOS
        % ================================================================
        
        rho_L5 = ...
            r_L5_all - r_FRT;
        
        d_L5 = ...
            vecnorm( ...
                rho_L5, ...
                2, ...
                2);
        
        u_L5 = ...
            rho_L5 ./ d_L5;
        
        H_L5 = [
            u_L5, ...
            ones(N_L5,1)
            ];

        %% ================================================================
        % Precompute individual H'H contributions
        % ================================================================

        M_NRHO = zeros(4,4,N_NRHO);
        M_L4   = zeros(4,4,N_L4);
        M_L5   = zeros(4,4,N_L5);


        for i = 1:N_NRHO

            h = H_NRHO(i,:)';

            M_NRHO(:,:,i) = ...
                h * h';

        end


        for i = 1:N_L4

            h = H_L4(i,:)';

            M_L4(:,:,i) = ...
                h * h';

        end


        for i = 1:N_L5

            h = H_L5(i,:)';

            M_L5(:,:,i) = ...
                h * h';

        end


        %% ================================================================
        % Evaluate all candidate combinations
        % ================================================================

        GDOP_candidates = ...
            Inf(nCombination,1);

        combination_index = 0;


        for NRHO_index = 1:N_NRHO

            M_EN = ...
                M_Earth + ...
                M_NRHO(:,:,NRHO_index);


            for L4_index = 1:N_L4

                M_EN4 = ...
                    M_EN + ...
                    M_L4(:,:,L4_index);


                for L5_index = 1:N_L5

                    combination_index = ...
                        combination_index + 1;


                    HTH = ...
                        M_EN4 + ...
                        M_L5(:,:,L5_index);


                    %% -------------------------------------------------
                    % Check geometry singularity
                    % -------------------------------------------------

                    if rcond(HTH) > GDOP_rcond_tolerance

                        Q = ...
                            HTH \ I4;

                        current_GDOP = ...
                            sqrt(trace(Q));

                    else

                        current_GDOP = Inf;

                    end


                    GDOP_candidates( ...
                        combination_index) = ...
                        current_GDOP;

                end
            end
        end


        %% ================================================================
        % Find Top-K candidates
        % ================================================================

        [sorted_GDOP, sorted_index] = ...
            mink(GDOP_candidates,TOP_K);


        %% -------------------------------------------------------------
        % Best solution
        % -------------------------------------------------------------

        best_index = ...
            sorted_index(1);

        best_GDOP_history( ...
            time_index,FRT_slot) = ...
            sorted_GDOP(1);

        best_combination_index_history( ...
            time_index,FRT_slot) = ...
            best_index;

        best_NRHO_ID_history( ...
            time_index,FRT_slot) = ...
            comb_NRHO(best_index);

        best_L4_ID_history( ...
            time_index,FRT_slot) = ...
            comb_L4(best_index);

        best_L5_ID_history( ...
            time_index,FRT_slot) = ...
            comb_L5(best_index);


        %% -------------------------------------------------------------
        % Top-K solutions
        % -------------------------------------------------------------

        for rank_index = 1:TOP_K

            candidate_index = ...
                sorted_index(rank_index);

            topK_GDOP_history( ...
                time_index, ...
                FRT_slot, ...
                rank_index) = ...
                sorted_GDOP(rank_index);

            topK_combination_index_history( ...
                time_index, ...
                FRT_slot, ...
                rank_index) = ...
                candidate_index;

            topK_NRHO_ID_history( ...
                time_index, ...
                FRT_slot, ...
                rank_index) = ...
                comb_NRHO(candidate_index);

            topK_L4_ID_history( ...
                time_index, ...
                FRT_slot, ...
                rank_index) = ...
                comb_L4(candidate_index);

            topK_L5_ID_history( ...
                time_index, ...
                FRT_slot, ...
                rank_index) = ...
                comb_L5(candidate_index);

        end

    end
end

%% ========================================================================
% Construct ID-aware FRT GDOP Lookup Table
% ========================================================================

nRows = ...
    nTime * N_FRT;


%% ------------------------------------------------------------------------
% Time / slot / ID information
% -------------------------------------------------------------------------

TimeIndex = ...
    repmat((1:nTime)',N_FRT,1);

Time = ...
    repmat(t_vector,N_FRT,1);

FRTSlot = ...
    repelem((1:N_FRT)',nTime);

FRT_ID = ...
    reshape( ...
        FRT_ID_history, ...
        [],1);


%% ------------------------------------------------------------------------
% FRT orbital phase
% -------------------------------------------------------------------------

PhaseAge = ...
    reshape( ...
        FRT_age_history, ...
        [],1);
T_FRT = tf_FRT;
NormalizedPhase = ...
    PhaseAge / T_FRT;


%% ------------------------------------------------------------------------
% FRT position
% -------------------------------------------------------------------------

X_matrix = ...
    squeeze(rv_FRT_dynamic(:,1,:));

Y_matrix = ...
    squeeze(rv_FRT_dynamic(:,2,:));

Z_matrix = ...
    squeeze(rv_FRT_dynamic(:,3,:));

X = reshape(X_matrix,[],1);
Y = reshape(Y_matrix,[],1);
Z = reshape(Z_matrix,[],1);


%% ------------------------------------------------------------------------
% Optimal GDOP information
% -------------------------------------------------------------------------

BestGDOP = ...
    reshape( ...
        best_GDOP_history, ...
        [],1);

CombinationIndex = ...
    reshape( ...
        best_combination_index_history, ...
        [],1);

BestNRHO_ID = ...
    reshape( ...
        best_NRHO_ID_history, ...
        [],1);

BestL4_ID = ...
    reshape( ...
        best_L4_ID_history, ...
        [],1);

BestL5_ID = ...
    reshape( ...
        best_L5_ID_history, ...
        [],1);


%% ========================================================================
% Final table
% ========================================================================

FRT_GDOP_lookup_table = table( ...
    TimeIndex, ...
    Time, ...
    FRTSlot, ...
    FRT_ID, ...
    PhaseAge, ...
    NormalizedPhase, ...
    X,Y,Z, ...
    BestGDOP, ...
    CombinationIndex, ...
    BestNRHO_ID, ...
    BestL4_ID, ...
    BestL5_ID);

end