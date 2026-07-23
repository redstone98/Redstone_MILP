%% RS_MILP_07_GS_to_SAT
function [x,z,x_hist,z_hist] = solve_ADMM(maxIters, rho, A_matrix, q, A_gj_info, b_gj_info, E, P_matrix, f_vector, gvec, E1_Gj, S_j_info, U_j_info)

%% Consensus ADMM with Gj-based local blocks

x = zeros(length(A_matrix(:,1)),1);
z = zeros(length(f_vector),1);

tol_primal = 1;
tol_dual   = 1;
verbose    = 1;

EP = E * P_matrix;

% ---- build Gj block index map ----
idx_block = cell(q,1);
for j = 1:q
    key = ['gs', num2str(j)];

    if isempty(E1_Gj.(key))
        idx_block{j} = [];
        continue;
    end

    U  = U_j_info.(key);
    Gj = S_j_info.(key);

    idx_block{j} = (U+1):(U+Gj);
end

% ---- dual variables ----
m = length(f_vector);
lambda = cell(q,1);
for j = 1:q
    lambda{j} = zeros(m,1);
end

x_hist = [];
z_hist = [];
r_hist = [];
s_hist = [];
cost_hist = [];

g = (gvec(:)).^2;

z_prev = z;

for k = 1:maxIters

    x_k = x;
    z_k = z;

    xbar_k = P_matrix' * x_k;

    z_gj_all = cell(q,1);
    xbar_new_all = cell(q,1);

    %% ========= (1) Local Gj updates
    for j = 1:q

        key = ['gs', num2str(j)];

        if isempty(E1_Gj.(key))
            z_gj_all{j} = z_k;
            xbar_new_all{j} = [];
            continue;
        end

        A_gj = A_gj_info.(key);
        b_gj = b_gj_info.(key);

        x_gj_k = E1_Gj.(key) * x_k;

        base0 = EP * xbar_k + f_vector - 1;

        U  = U_j_info.(key);
        Gj = S_j_info.(key);

        EP_j = EP(:, (U+1):(U+Gj));

        params = struct();
        params.OutputFlag = 0;

        sol = solve_local_admm_milp_gurobi( ...
            A_gj, b_gj, EP_j, base0, x_gj_k, ...
            g, z_k, lambda{j}, rho, params);

        if ~sol.success
            if verbose
                fprintf('[ADMM-Gj] iter=%d gs%d local solve failed (status=%s). Keep.\n', ...
                    k, j, sol.status);
            end

            z_gj_all{j} = z_k;
            xbar_new_all{j} = x_gj_k;
            continue;
        end

        x_gj_new = sol.x_gj;   % solver 내부 변수명이 x_si여도 값은 Gj block
        z_gj_new = sol.z_gj;

        z_gj_all{j} = z_gj_new;
        xbar_new_all{j} = x_gj_new;
    end

    %% ========= (2) Global x update
    xbar_seq = xbar_k;

    for j = 1:q
        if isempty(idx_block{j})
            continue;
        end

        key = ['gs', num2str(j)];

        U  = U_j_info.(key);
        Gj = S_j_info.(key);

        xbar_seq(U+1:U+Gj) = xbar_new_all{j};
    end

    x = P_matrix * xbar_seq;

    %% ========= (3) Global z update
    active = false(q,1);
    for j = 1:q
        key = ['gs', num2str(j)];
        active(j) = ~isempty(E1_Gj.(key));
    end

    Aset = find(active);
    cnt  = numel(Aset);

    if cnt == 0
        z = z_k;
    else
        xbar = P_matrix' * x;

        lb1 = EP * xbar + f_vector - 1;
        lb2 = zeros(m,1);
        lb  = max(lb1, lb2);

        v = zeros(m,1);
        for jj = 1:cnt
            j = Aset(jj);
            v = v + z_gj_all{j} + (1/rho) * lambda{j};
        end
        v = v / cnt;

        z = max(v, lb);
    end

    %% ========= (4) Dual update
    for j = 1:q
        key = ['gs', num2str(j)];

        if isempty(E1_Gj.(key))
            continue;
        end

        lambda{j} = lambda{j} + rho * (z_gj_all{j} - z);
    end

    %% ========= (5) Residual checks
    r_norm = 0;

    for j = 1:q
        key = ['gs', num2str(j)];

        if isempty(E1_Gj.(key))
            continue;
        end

        r_norm = r_norm + norm(z_gj_all{j} - z, 2)^2;
    end

    r_norm = sqrt(r_norm);

    s_norm = norm(z - z_prev, 2);

    cost_hist = [cost_hist, g' * z];

    x_hist = [x_hist, x];
    z_hist = [z_hist, z];
    r_hist = [r_hist, r_norm];
    s_hist = [s_hist, s_norm];

    if verbose
        fprintf('[ADMM-Gj] iter=%3d cost = %.3e |r|=%.3e |s|=%.3e\n', ...
            k, g' * z, r_norm, s_norm);
    end

    if (r_norm <= tol_primal) && (s_norm <= tol_dual)
        if verbose
            fprintf('[ADMM-Gj] Converged at iter=%d\n', k);
        end
        break;
    end

    z_prev = z;
end

end