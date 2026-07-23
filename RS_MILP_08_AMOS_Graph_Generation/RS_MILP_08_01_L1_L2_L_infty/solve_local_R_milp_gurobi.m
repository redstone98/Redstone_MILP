function sol = solve_local_R_milp_gurobi(A_si, b_si, EP_i, base0, x_si_k, gvec, params)
% Local BCD MILP:
%
%   min_{x_si, R_si} R_si
%
%   s.t. A_si x_si <= b_si
%
%        R_si * 1 >= diag(gvec) *
%        (base0 + EP_i*(x_si - x_si_k))
%
% where
%   base0 = EP*xbar_k + f_vector - 1
    sol = struct('success', false, 'x_si', [], 'R', [], 'status', '');

    if nargin < 7 || isempty(params)
        params = struct();
    end

    % -----------------------------
    % Preprocessing
    % -----------------------------
    A_si   = sparse(A_si);
    EP_i   = sparse(EP_i);
    nSi = size(A_si,2);
    m   = length(base0);

    if length(gvec) ~= m
        error('length(gvec) must match length(base0).');
    end

    % -----------------------------
    % Variables: v = [x_si; R_si]
    % -----------------------------
    nvar = nSi + 1;

    model.modelsense = 'min';
    model.obj = [zeros(nSi,1); 1];

    model.vtype = [repmat('B', nSi, 1); 'C'];
    model.lb = [zeros(nSi,1); 0];
    model.ub = [ones(nSi,1); inf];

    % ======================================================
    % Constraint 1:
    % A_si x_si <= b_si
    % ======================================================
    A1 = sparse(size(A_si,1), nvar);
    A1(:,1:nSi) = A_si;

    rhs1 = b_si;
    sense1 = repmat('<', size(A_si,1), 1);

    % ======================================================
    % Constraint 2:
    %
    % R_si*1 >= diag(gvec)*(base0 + EP_i*(x_si - x_si_k))
    %
    % Expand:
    %
    % R_si*1 >= diag(gvec)*(base0 - EP_i*x_si_k + EP_i*x_si)
    %
    % Equivalent Gurobi <= form:
    %
    % diag(gvec)*EP_i*x_si - R_si*1
    % <= -diag(gvec)*(base0 - EP_i*x_si_k)
    % ======================================================
    base_i = base0 - EP_i*x_si_k;

    A2 = sparse(m, nvar);
    A2(:,1:nSi) = spdiags(gvec, 0, m, m) * EP_i;
    A2(:,end)   = -ones(m,1);

    rhs2 = -gvec .* base_i;
    sense2 = repmat('<', m, 1);

    % -----------------------------
    % Build Gurobi model
    % -----------------------------
    model.A = sparse([A1; A2]);
    model.rhs = [rhs1; rhs2];
    model.sense = [sense1; sense2];

    % -----------------------------
    % Solve
    % -----------------------------
    try
        result = gurobi(model, params);
    catch ME
        sol.status = ['gurobi_error: ', ME.message];
        return;
    end

    sol.status = result.status;

    if strcmpi(result.status,'OPTIMAL') || strcmpi(result.status,'SUBOPTIMAL')
        v = result.x;
        sol.x_si = full(double(v(1:nSi)));
        sol.R    = full(double(v(end)));
        sol.success = true;
    end
end