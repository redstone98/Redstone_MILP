function [x,z, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, p, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info)

%% =========================
%  BCD for MILP with z >= E*P*xbar + f - 1
%  - each agent i updates x_{s_i} (binary) sequentially
%  - z is continuous and is re-optimized in each local subproblem
%% =========================
N = length(A_matrix(:,1));


% ---- user provided initializations (keep) ----
x = zeros(length(A_matrix(:,1)),1);
z = zeros(length(f_vector),1);

x_history_BCD = [];
z_history_BCD = [];

% ---- BCD options ----
maxOuterIters = 50;          % outer sweeps (full i=1..p)
verbose       = 1;

% ---- precompute constant matrices ----
% P_matrix: permutation such that x = P * [x_s1; ...; x_sp]
% EP = E * P  (used in all agents)
EP = E * P_matrix;

% ---- build block index map (idx_i for each agent) ----
% You already have:
%   U = sum_{j<i} |S_j|,  S_i is list (or just length defines block size)
% We'll build idx_i = (U+1):(U+|S_i|)
idx_block = cell(p,1);
for i = 1:p
    if isempty(E1_Si.(['sat', num2str(i)])) 
        idx_block{i} = [];
        continue;
    end
    U_i = U_i_info.(['sat', num2str(i)]);
    S_i = S_i_info.(['sat', num2str(i)]);
    nSi = length(S_i);
    idx_block{i} = (U_i+1):(U_i+nSi);
end

BCD_handle = 1;
outerIter  = 0;

while BCD_handle == 1
    outerIter = outerIter + 1;
    if outerIter > maxOuterIters
        if verbose, fprintf('[BCD] Reached maxOuterIters=%d\n', maxOuterIters); end
        break;
    end

    x_old_sweep = x;
    z_old_sweep = z;

    for i = 1:p
        % skip empty agent
        if isempty(E1_Si.(['sat', num2str(i)]))
            continue;
        end

        % ----- snapshot current iterate -----
        x_k = x;
        z_k = z; %#ok<NASGU>  % (kept if you want warm-start logic)

        % ----- local constraints -----
        A_si = A_si_info.(['sat', num2str(i)]);
        b_si = b_si_info.(['sat', num2str(i)]);

        % ----- stacked vector at k: xbar^k = P' * x^k -----
        xbar_k = P_matrix' * x_k;

        % ----- identify current local block -----
        x_si_k = E1_Si.(['sat', num2str(i)]) * x_k;       % this is the "current" local block

        % ----- coupling base term: base0 = EP*xbar_k + f - 1 -----
        base0  = EP * xbar_k + f_vector - 1;

        % ----- the only part that multiplies x_si variable -----

        U_i = U_i_info.(['sat', num2str(i)]);
        S_i = S_i_info.(['sat', num2str(i)]);
        EP_i = EP(:, (U_i+1):(U_i+S_i));          % == EP * J_i

        J_i = sparse(N, S_i);              % (U+|S_i|+V) x |S_i|
        J_i(U_i+1:U_i+S_i, :) = speye(S_i);    % middle block = I_|S_i|

        % ----- Build and solve local MIQP/MILP -----
        % Variables: [x_si; z]
        %   x_si in {0,1}^{|S_i|},  z in R^{m}, typically with z >= 0
        %
        % Constraints:
        %   (1) A_si * x_si <= b_si
        %   (2) z >= base0 + EP_i*(x_si - x_si_k)
        %       -> -z + EP_i*x_si <= -(base0 - EP_i*x_si_k)
        %   (3) z >= 0   (recommended if z is slack-like)
        %
        % Objective:
        %   min sum_j G_j * z_j^2   (i.e., z' diag(G) z)
        %
        % NOTE: set Gvec (length m) as your diag weights.
        % If you have G as matrix, pass diag(G).

        % choose your weight vector for quadratic objective:
        % - if you already have G as vector: use it directly
        % - if you have G as diagonal matrix: use diag(G)
        g = gvec.^2;


         % solve local subproblem with MIQP
        params = optimoptions('intlinprog','Display','off');

        sol = solve_local_milp_intlinprog(A_si, b_si, EP_i, base0, x_si_k, g, params);

        % if infeasible or failed, skip update gracefully
        if ~sol.success
            if verbose
                fprintf('[BCD] sat%d local solve failed (status=%s). Keep x,z.\n', i, sol.status);
            end
            continue;
        end

        x_si_new = sol.x_si;
        z_new    = sol.z;

        % ----- write back to global x via stacked vector -----
        xbar_new = xbar_k + J_i* (x_si_new - x_si_k);
        x = P_matrix * xbar_new;
        z = z_new;

        % history
        x_history_BCD = [x_history_BCD, x]; %#ok<AGROW>
        z_history_BCD = [z_history_BCD, z]; %#ok<AGROW>
    end

    % ----- stopping check (end of sweep) -----
    if isequal(x, x_old_sweep) && norm(z - z_old_sweep, Inf) <= 1e-12
        if verbose
            fprintf('[BCD] Converged at outerIter=%d (no change after full sweep)\n', outerIter);
        end
        break;
    else
        if verbose
            dx = nnz(x ~= x_old_sweep);
            dz = norm(z - z_old_sweep, Inf);
            fprintf('[BCD] outerIter=%d done. nnz(dx)=%d, ||dz||_inf=%.3e\n', outerIter, dx, dz);
        end
    end
end

end