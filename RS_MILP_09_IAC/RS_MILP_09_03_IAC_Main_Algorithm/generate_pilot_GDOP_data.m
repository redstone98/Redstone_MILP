function [rv_FRT_one, rv_NRHO_constellation, rv_L4_constellation, rv_L5_constellation, nCombination, GDOP_history, NRHO_index_history, L4_index_history, L5_index_history, combination_table] = ...
    generate_pilot_GDOP_data(t_vector, rv0_FRT_SATs,rv0_NRHO_south_SATs, rv0_L4_vertical_SATs, rv0_L5_vertical_SATs)

global mu

number_of_NRHO_SATs = length(rv0_NRHO_south_SATs(:,1));
number_of_L4_SATs = length(rv0_L4_vertical_SATs(:,1));
number_of_L5_SATs = length(rv0_L5_vertical_SATs(:,1));


rv0_FRT = rv0_FRT_SATs(1,:)';

odeopt = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t, rv_FRT_one] = ode45(@CRTBP, t_vector, rv0_FRT, odeopt);

rv_NRHO_constellation = zeros(length(t_vector),6,number_of_NRHO_SATs);
rv_L4_constellation =  zeros(length(t_vector),6,number_of_L4_SATs);
rv_L5_constellation = zeros(length(t_vector),6,number_of_L5_SATs);

for NRHO_index = 1:number_of_NRHO_SATs
  rv0_NRHO_i = rv0_NRHO_south_SATs(NRHO_index,:)';
  [t, rv_NRHO_i] = ode45(@CRTBP, t_vector, rv0_NRHO_i, odeopt);  
  rv_NRHO_constellation(:,:,NRHO_index) = rv_NRHO_i;
end


for L4_index = 1:number_of_L4_SATs
  rv0_L4_i = rv0_L4_vertical_SATs(L4_index,:)';
  [t, rv_L4_i] = ode45(@CRTBP, t_vector, rv0_L4_i, odeopt);  
  rv_L4_constellation(:,:,L4_index) = rv_L4_i;
end


for L5_index = 1:number_of_L5_SATs
  rv0_L5_i = rv0_L5_vertical_SATs(L5_index,:)';
  [t, rv_L5_i] = ode45(@CRTBP, t_vector, rv0_L5_i, odeopt);  
  rv_L5_constellation(:,:,L5_index) = rv_L5_i;
end

%% Calculate GDOP
%
% Navigation sources:
%   1) Earth
%   2) One NRHO satellite
%   3) One L4 satellite
%   4) One L5 satellite
%
% For each timestep:
% number of combinations =
% number_of_NRHO_SATs * number_of_L4_SATs * number_of_L5_SATs

nTime = length(t_vector);

nCombination = ...
    number_of_NRHO_SATs * ...
    number_of_L4_SATs * ...
    number_of_L5_SATs;

%% ------------------------------------------------------------------------
% Preallocation
% -------------------------------------------------------------------------

GDOP_history = NaN(nTime, nCombination);

% Satellite index tracking
NRHO_index_history = zeros(nCombination,1);
L4_index_history   = zeros(nCombination,1);
L5_index_history   = zeros(nCombination,1);

% Combination index
combination_index = 0;

for NRHO_index = 1:number_of_NRHO_SATs
    for L4_index = 1:number_of_L4_SATs
        for L5_index = 1:number_of_L5_SATs

            combination_index = combination_index + 1;

            NRHO_index_history(combination_index) = NRHO_index;
            L4_index_history(combination_index)   = L4_index;
            L5_index_history(combination_index)   = L5_index;

        end
    end
end

%% ------------------------------------------------------------------------
% GDOP calculation
% -------------------------------------------------------------------------

% Earth position
r_Earth = [-mu; 0; 0];

for time_index = 1:nTime

    % -------------------------------------------------------------
    % FRT client satellite position
    % -------------------------------------------------------------
    r_FRT = rv_FRT_one(time_index,1:3)';

    combination_index = 0;

    for NRHO_index = 1:number_of_NRHO_SATs

        r_NRHO = ...
            rv_NRHO_constellation(time_index,1:3,NRHO_index)';

        for L4_index = 1:number_of_L4_SATs

            r_L4 = ...
                rv_L4_constellation(time_index,1:3,L4_index)';

            for L5_index = 1:number_of_L5_SATs

                combination_index = combination_index + 1;

                r_L5 = ...
                    rv_L5_constellation(time_index,1:3,L5_index)';

                %% -------------------------------------------------
                % Line-of-sight vectors
                % -------------------------------------------------

                rho_Earth = r_Earth - r_FRT;
                rho_NRHO  = r_NRHO  - r_FRT;
                rho_L4    = r_L4    - r_FRT;
                rho_L5    = r_L5    - r_FRT;

                %% Unit LOS vectors
                u_Earth = rho_Earth / norm(rho_Earth);
                u_NRHO  = rho_NRHO  / norm(rho_NRHO);
                u_L4    = rho_L4    / norm(rho_L4);
                u_L5    = rho_L5    / norm(rho_L5);

                %% -------------------------------------------------
                % Geometry Matrix
                %
                % Position + clock bias estimation
                % -------------------------------------------------

                H = [
                    u_Earth'  1
                    u_NRHO'   1
                    u_L4'     1
                    u_L5'     1
                ];

                %% -------------------------------------------------
                % GDOP
                % -------------------------------------------------

                HTH = H' * H;

                % Numerical singularity check
                if rcond(HTH) > 1e-12

                    Q = inv(HTH);

                    GDOP = sqrt(trace(Q));

                else

                    % Bad / singular geometry
                    GDOP = Inf;

                end

                GDOP_history(time_index, combination_index) = GDOP;

            end
        end
    end
end

combination_table = table( ...
    (1:nCombination)', ...
    NRHO_index_history, ...
    L4_index_history, ...
    L5_index_history, ...
    'VariableNames', { ...
        'CombinationIndex', ...
        'NRHO_Index', ...
        'L4_Index', ...
        'L5_Index'})

end
