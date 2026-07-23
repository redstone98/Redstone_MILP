%% RS_MILP_07_GS_to_SAT
function sol = solve_local_admm_milp_gurobi(A_gj, b_gj, EP_j, base0, x_gj_k, gvec, z_k, lambda_j, rho, params)
% Solve local Gj-based MILP/MIQP:
%
%   min_{x_gj in {0,1}^n, z in R^m}
%       g' z + lambda_j' (z - z_k) + (rho/2) ||z - z_k||_2^2
%
%   s.t.
%       A_gj x_gj <= b_gj
%       z >= base0 + EP_j*(x_gj - x_gj_k)
%       z >= 0
%
% Variables v = [x_gj; z].

    if nargin < 10 || isempty(params)
        params = struct();
    end

    if ~isfield(params,'OutputFlag')
        params.OutputFlag = 0;
    end

    if ~isfield(params,'LogToConsole')
        params.LogToConsole = 0;
    end

    if ~isfield(params,'LogFile')
        params.LogFile = '';
    end

    gvec     = gvec(:);
    z_k      = z_k(:);
    lambda_j = lambda_j(:);
    base0    = base0(:);
    x_gj_k   = x_gj_k(:);

    n = size(A_gj,2);      % |x_gj|
    m = length(base0);     % |z|

    % -------------------------
    % Objective:
    %
    % g'z + lambda_j'(z - z_k) + (rho/2)||z - z_k||^2
    %
    % = (rho/2) z'z + (g + lambda_j - rho z_k)' z + const
    %
    % Gurobi form:
    %   min (1/2) v'Qv + c'v
    % -------------------------
    Q = sparse(n+m, n+m);
    Q(n+1:n+m, n+1:n+m) = rho * speye(m);

    c = zeros(n+m,1);
    c(n+1:n+m) = gvec + lambda_j - rho*z_k;

    % -------------------------
    % Constraints:
    %
    % (1) A_gj x_gj <= b_gj
    %
    % (2) z >= base0 + EP_j*(x_gj - x_gj_k)
    %
    %     z >= base0 + EP_j*x_gj - EP_j*x_gj_k
    %
    %     EP_j*x_gj - z <= -base0 + EP_j*x_gj_k
    %
    %     [EP_j, -I] [x_gj; z] <= -(base0 - EP_j*x_gj_k)
    %
    % (3) z >= 0 handled by lower bound
    % -------------------------
    A1   = [A_gj, sparse(size(A_gj,1), m)];
    rhs1 = b_gj(:);

    A2   = [EP_j, -speye(m)];
    rhs2 = -(base0 - EP_j*x_gj_k);

    A    = sparse([A1; A2]);
    rhs  = [rhs1; rhs2];

    sense = repmat('<', size(A,1), 1);

    % -------------------------
    % Bounds
    % -------------------------
    lb = -inf(n+m,1);
    ub =  inf(n+m,1);

    lb(1:n) = 0;
    ub(1:n) = 1;

    lb(n+1:n+m) = 0;

    % -------------------------
    % Variable types
    % -------------------------
    vtype = repmat('C', n+m, 1);
    vtype(1:n) = 'B';

    % -------------------------
    % Gurobi model
    % -------------------------
    model.Q = Q;
    model.obj = c;
    model.A = A;
    model.rhs = rhs;
    model.sense = sense;
    model.lb = lb;
    model.ub = ub;
    model.vtype = vtype;
    model.modelsense = 'min';

    try
        result = gurobi(model, params);
    catch ME
        sol.success = false;
        sol.status  = ['error: ' ME.message];
        sol.x_gj = [];
        sol.z_gj = [];
        return;
    end

    sol.status = result.status;

    if isfield(result,'x') && strcmpi(result.status,'OPTIMAL')
        v = result.x;

        sol.x_gj = v(1:n);
        sol.z_gj = v(n+1:n+m);
        sol.success = true;
    else
        sol.success = false;
        sol.x_gj = [];
        sol.z_gj = [];
    end
end