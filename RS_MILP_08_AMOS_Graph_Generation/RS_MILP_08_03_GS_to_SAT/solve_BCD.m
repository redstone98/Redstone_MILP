%% RS_MILP_07_GS_to_SAT
function [x,z, x_history_BCD, z_history_BCD] = solve_BCD(A_matrix, q, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec ,E1_Gj, S_j_info, U_j_info)

%% =========================
%  BCD for MILP with z >= E*P*xbar + f - 1
%  - each ground-point / ground-station block j updates x_{g_j} sequentially
%  - x_{g_j} is the binary contact-selection vector associated with Gj block
%  - z is continuous and is re-optimized in each local subproblem
%% =========================
N = length(A_matrix(:,1));


% ---- user provided initializations ----
x = zeros(N,1);
z = zeros(length(f_vector),1);

x_history_BCD = [];
z_history_BCD = [];

% ---- BCD options ----
maxOuterIters = 50;          % outer sweeps over all ground blocks j = 1,...,q
verbose       = 1;

% ---- precompute constant matrices ----
% P_matrix: permutation such that
%   x = P_matrix * [x_g1; x_g2; ...; x_gq]
%
% EP = E * P_matrix
% This is used to express the global coupling term in the Gj-stacked order.
EP = E * P_matrix;

% ---- build block index map for each Gj block ----
% U_j is the starting offset of block j in the stacked vector
%   xbar = [x_g1; x_g2; ...; x_gq]
%
% S_j_info.(['gs',num2str(j)]) gives the size or index list of Gj block.
% idx_block{j} identifies the position of x_gj inside xbar.
idx_block = cell(q,1);
for j = 1:q
    if isempty(E1_Gj.(['gs', num2str(j)])) 
        idx_block{j} = [];
        continue;
    end
    U_j = U_j_info.(['gs', num2str(j)]);
    S_j = S_j_info.(['gs', num2str(j)]);
    nSj = length(S_j);
    idx_block{j} = (U_j+1):(U_j+nSj);
end

BCD_handle = 1;
outerIter  = 0;

while BCD_handle == 1
    outerIter = outerIter + 1;
    if outerIter > maxOuterIters
        if verbose, fprintf('[BCD] Reached maxOuterIters=%d\n', maxOuterIters); end
        break;
    end

    z_old_sweep = z;

    for j = 1:q
        % skip empty Gj block
        if isempty(E1_Gj.(['gs', num2str(j)]))
            continue;
        end

        % ----- snapshot current iterate -----
        x_k = x;
        z_k = z; %#ok<NASGU>  % kept for possible warm-start logic

        % ----- local Gj constraints -----
        % A_gj x_gj <= b_gj
        A_gj = A_gj_info.(['gs', num2str(j)]);
        b_gj = b_gj_info.(['gs', num2str(j)]);

        % ----- current stacked vector -----
        % xbar_k is the Gj-stacked version of the current global vector x_k:
        %
        %   xbar_k = [x_g1^k; x_g2^k; ...; x_gq^k]
        %
        xbar_k = P_matrix' * x_k;

        % ----- current local Gj block -----
        % Extract x_gj^k from the global decision vector x_k.
        x_gj_k = E1_Gj.(['gs', num2str(j)]) * x_k;

        % ----- coupling base term -----
        % For current xbar_k:
        %
        %   base0 = EP*xbar_k + f_vector - 1
        %
        % During the local update, only x_gj is allowed to change.
        base0  = EP * xbar_k + f_vector - 1;

        % ----- local coupling matrix for Gj block -----
        % EP_j contains only the columns of EP corresponding to x_gj.
        %
        % Therefore:
        %
        %   EP*xbar_new + f - 1
        %   = base0 + EP_j*(x_gj_new - x_gj_k)
        %
        U_j = U_j_info.(['gs', num2str(j)]);
        S_j = S_j_info.(['gs', num2str(j)]);
        EP_j = EP(:, (U_j+1):(U_j+S_j));

        % ----- injection matrix for writing x_gj back into xbar -----
        % J_i maps a local Gj update into the stacked global vector xbar.
        %
        %   xbar_new = xbar_k + J_i*(x_gj_new - x_gj_k)
        %
        % Although named J_i here, it represents the injection matrix for
        % the current Gj block.
        J_i = sparse(N, S_j);
        J_i(U_j+1:U_j+S_j, :) = speye(S_j);

        % ----- Build and solve local Gj MILP/MIQP -----
        % Variables:
        %   [x_gj; z]
        %
        % where:
        %   x_gj in {0,1}^{|G_j|}
        %   z    in R^m
        %
        % Constraints:
        %   (1) A_gj * x_gj <= b_gj
        %
        %   (2) z >= base0 + EP_j*(x_gj - x_gj_k)
        %
        %       equivalently:
        %
        %       EP_j*x_gj - z <= -(base0 - EP_j*x_gj_k)
        %
        %   (3) z >= 0
        %
        % Objective:
        %   min z' diag(g) z
        %
        % where:
        %   g = gvec.^2
        %
        g = gvec.^2;

        % solve local subproblem
        params = optimoptions('intlinprog','Display','off');

        sol = solve_local_milp_intlinprog(A_gj, b_gj, EP_j, base0, x_gj_k, g, params);

        % if infeasible or failed, skip update gracefully
        if ~sol.success
            if verbose
                fprintf('[BCD] GS%d local solve failed (status=%s). Keep x,z.\n', j, sol.status);
            end
            continue;
        end

        x_gj_new = sol.x_gj;
        z_new    = sol.z;

        % ----- write updated Gj block back to global x -----
        % First update the stacked vector xbar, then map it back to the
        % original global ordering:
        %
        %   xbar_new = xbar_k + J_i*(x_gj_new - x_gj_k)
        %   x        = P_matrix*xbar_new
        %
        xbar_new = xbar_k + J_i* (x_gj_new - x_gj_k);
        x = P_matrix * xbar_new;
        z = z_new;

        % history
        x_history_BCD = [x_history_BCD, x]; %#ok<AGROW>
        z_history_BCD = [z_history_BCD, z]; %#ok<AGROW>
    end

    % ----- stopping check at the end of one full Gj sweep -----
    if norm(z - z_old_sweep, Inf) <= 1e-12
        if verbose
            fprintf('[BCD] Converged at outerIter=%d (no change after full sweep)\n', outerIter);
        end
        break;
    else
        if verbose
            dz = norm(z - z_old_sweep, Inf);
            fprintf('[BCD] outerIter=%d done.  ||dz||_inf=%.3e\n', outerIter, dz);
        end
    end
end

end