function [x,z, x_hist, z_hist] = solve_ADMM(maxIters, rho, A_matrix, p, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info)


%% =========================
%  Consensus ADMM for (MIQP/MILP-like) with z >= E*P*xbar + f - 1
%   - each agent i updates (x_{s_i}, z_{s_i}) in parallel/serial
%   - global z enforces consensus: z_{s_i} = z
%   - dual variables lambda_i (same dimension as z)
%% =========================

% ---- user provided initializations (keep) ----
x = zeros(length(A_matrix(:,1)),1);     % global x (binary)
z = zeros(length(f_vector),1);          % global consensus z (continuous)

% ---- ADMM options ----
tol_primal = 1;      % ||r|| tolerance
tol_dual   = 1;      % ||s|| tolerance
verbose    = 1;

% ---- precompute constant matrices ----
EP = E * P_matrix;      % used for z lower bound

% ---- build block index map (idx_i for each agent) ----
idx_block = cell(p,1);
for i = 1:p
    if isempty(E1_Si.(['sat', num2str(i)]))
        idx_block{i} = [];
        continue;
    end
    U = U_i_info.(['sat', num2str(i)]);
    S_i = S_i_info.(['sat', num2str(i)]);
    nSi = S_i;
    idx_block{i} = (U+1):(U+nSi);
end

% ---- dual variables (one lambda per agent, dimension = length(z)) ----
m = length(f_vector);
lambda = cell(p,1);
for i = 1:p
    lambda{i} = zeros(m,1);
end

% ---- histories (optional) ----
x_hist = [];
z_hist = [];
r_hist = [];
s_hist = [];

cost_hist = [];

% ---- weight for quadratic on z_{s_i} ----
g = (gvec(:)).^2;   % objective uses z' diag(g) z

% ---- main ADMM loop ----
z_prev = z;

for k = 1:maxIters

    % snapshot
    x_k = x;
    z_k = z;

    % stacked xbar = P' x
    xbar_k = P_matrix' * x_k;

    % store local solutions
    z_si_all = cell(p,1);
    xbar_new_all = cell(p,1);

    % ========= (1) Local updates: for each agent i, solve MIQP in (x_si, z_si)
    for i = 1:p
        if isempty(E1_Si.(['sat', num2str(i)]))
            z_si_all{i} = z_k;          % dummy
            xbar_new_all{i} = xbar_k;   % dummy
            continue;
        end

        % local constraints
        A_si = A_si_info.(['sat', num2str(i)]);
        b_si = b_si_info.(['sat', num2str(i)]);

        % current local block (from global x)
        x_si_k = E1_Si.(['sat', num2str(i)]) * x_k;

        % coupling base term: base0 = EP*xbar_k + f - 1
        base0 = EP * xbar_k + f_vector - 1;

        % EP_i: columns corresponding to agent i block inside xbar
        U  = U_i_info.(['sat', num2str(i)]);
        S  = S_i_info.(['sat', num2str(i)]);
        EP_i = EP(:, (U+1):(U+S));   % m x |S_i|

        % solve local MIQP:
        %   min_{x_si in {0,1}^n, z in R^m}
        %       g' z + lambda' (z - z_k) + (rho/2) ||z - z_k||_2^2
        %   s.t. A_si x_si <= b_si, x_si binary
        %        z_si >= base0 + EP_i*(x_si - x_si_k)
        %        z_si >= 0
        params = struct();
        params.OutputFlag = 0;

        sol = solve_local_admm_milp_gurobi(A_si, b_si, EP_i, base0, x_si_k, g, z_k, lambda{i}, rho, params);

        if ~sol.success
            if verbose
                fprintf('[ADMM] iter=%d sat%d local solve failed (status=%s). Keep.\n', k, i, sol.status);
            end
            z_si_all{i} = z_k;
            xbar_new_all{i} = xbar_k;
            continue;
        end

        x_si_new = sol.x_si;
        z_si_new = sol.z_si;
        z_si_all{i} = z_si_new;
        xbar_new_all{i} = x_si_new;
    end

    % ========= (2) Global x update
    % 보통 consensus ADMM에서는 x를 "각 agent의 x_si를 합쳐" 업데이트해야 합니다.
    % 여기서는 가장 단순히 "마지막 sweep의 xbar_new_all{i}"를 사용하거나,
    % 혹은 다수결/평균 후 rounding 등을 설계할 수 있습니다.
    %
    % 아래는 "sequential sweep"처럼 i=1..p 순서로 반영한 결과를 사용:
    xbar_seq = xbar_k;
    for i = 1:p
        if isempty(idx_block{i}), continue; end
        U  = U_i_info.(['sat', num2str(i)]);
        S  = S_i_info.(['sat', num2str(i)]);
        xbar_seq(U+1:U+S) = xbar_new_all{i};
    end
    x = P_matrix * xbar_seq;

    % ========= (3) Global z update (Option ii)
    m = length(f_vector);
    
    % active agent set
    active = false(p,1);
    for i = 1:p
        active(i) = ~isempty(E1_Si.(['sat', num2str(i)]));
    end
    Aset = find(active);
    cnt  = numel(Aset);
    
    if cnt == 0
        z = z_k;
    else
        % compute coupling lower bound using UPDATED x (after step (2))
        xbar = P_matrix' * x;                 
        lb1  = EP * xbar + f_vector - 1;
        lb2  = zeros(m,1);
        lb   = max(lb1, lb2);                 % element-wise lower bound
    
        % build "center" v = (1/cnt) * sum_i ( z_si + (1/rho)*lambda_i )
        v = zeros(m,1);
        for ii = 1:cnt
            i = Aset(ii);
            v = v + ( z_si_all{i} + (1/rho) * lambda{i} );
        end
        v = v / cnt;
    
        % projection onto z >= lb  (element-wise)
        z = max(v, lb);
    end



    % ========= (4) Dual update
    for i = 1:p
        if isempty(E1_Si.(['sat', num2str(i)])), continue; end
        lambda{i} = lambda{i} + rho*( z_si_all{i} - z );
    end

    % ========= (5) Residual checks
    % primal residual: r_i = z_si - z
    r_norm = 0;
    cnt = 0;
    for i = 1:p
        if isempty(E1_Si.(['sat', num2str(i)])), continue; end
        r_norm = r_norm + norm(z_si_all{i} - z, 2)^2;
        cnt = cnt + 1;
    end
    r_norm = sqrt(r_norm);

    % dual residual: s = rho*(z - z_prev)
    s_norm = norm(z - z_prev, 2);


    % Calculate cost
    cost_hist = [cost_hist, g'*z];

    % histories
    x_hist = [x_hist, x];
    z_hist = [z_hist, z];
    r_hist = [r_hist, r_norm];
    s_hist = [s_hist, s_norm];


    if verbose
        fprintf('[ADMM] iter=%3d cost = %.3e |r|=%.3e |s|=%.3e\n', k, g'*z, r_norm, s_norm);
    end

    if (r_norm <= tol_primal) && (s_norm <= tol_dual)
        if verbose, fprintf('[ADMM] Converged at iter=%d\n', k); end
        break;
    end

    z_prev = z;
end


end