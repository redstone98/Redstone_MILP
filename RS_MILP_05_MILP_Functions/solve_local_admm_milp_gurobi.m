function sol = solve_local_admm_milp_gurobi(A_si, b_si, EP_i, base0, x_si_k, gvec, z_k, lambda_i, rho, params)
% Solve local MILP:
%
%   min_{x_si in {0,1}^n, z in R^m}
%       g' z + lambda' (z - z_k) + (rho/2) ||z - z_k||_2^2
%   s.t.
%       A_si x_si <= b_si
%       z >= base0 + EP_i*(x_si - x_si_k)
%       z >= 0
%
% Variables v = [x; z].

    if nargin < 10 || isempty(params), params = struct(); end
    if ~isfield(params,'OutputFlag'),   params.OutputFlag = 0; end
    if ~isfield(params,'LogToConsole'), params.LogToConsole = 0; end
    if ~isfield(params,'LogFile'),      params.LogFile = ''; end

    gvec     = gvec(:);
    z_k      = z_k(:);
    lambda_i = lambda_i(:);

    n = size(A_si,2);      % |x_si|
    m = length(base0);     % |z|

    % -------------------------
    % Objective expansion:
    % g'z + lambda'(z - z_k) + (rho/2)||z - z_k||^2
    %
    % = (rho/2) z'z + (g + lambda - rho z_k)' z + const
    %
    % Gurobi form: min (1/2) v'Qv + c'v
    % -> put Qzz = rho*I  (because (1/2) z'(rho I) z = (rho/2) z'z)
    % -> linear term on z: c_z = g + lambda - rho z_k
    % -------------------------
    Q = sparse(n+m, n+m);
    Q(n+1:n+m, n+1:n+m) = rho * speye(m);

    c = zeros(n+m,1);
    c(n+1:n+m) = gvec + lambda_i - rho*z_k;

    % -------------------------
    % Constraints:
    % (1) A_si x <= b_si
    % (2) z >= base0 + EP_i*(x - x_k)
    %     -> -z + EP_i*x <= -(base0 - EP_i*x_k)
    % (3) z >= 0 handled by lb on z
    % -------------------------
    A1   = [A_si, sparse(size(A_si,1), m)];
    rhs1 = b_si;

    A2   = [EP_i, -speye(m)];
    rhs2 = -(base0 - EP_i*x_si_k);

    A    = sparse([A1; A2]);
    rhs  = [rhs1; rhs2];
    sense = repmat('<', size(A,1), 1);

    % bounds
    lb = -inf(n+m,1);
    ub =  inf(n+m,1);

    lb(1:n) = 0;   ub(1:n) = 1;     % x binary bounds
    lb(n+1:n+m) = 0;                % z >= 0

    vtype = repmat('C', n+m, 1);
    vtype(1:n) = 'B';

    model.Q = Q;
    model.obj = c;
    model.A = A;
    model.rhs = rhs;
    model.sense = sense;
    model.lb = lb;
    model.ub = ub;
    model.vtype = vtype;
    model.modelsense = 'min';

    % Solve
    try
        result = gurobi(model, params);
    catch ME
        sol.success = false;
        sol.status  = ['error: ' ME.message];
        sol.x_si = [];
        sol.z_si = [];
        return;
    end

    sol.status = result.status;

    if isfield(result,'x') && strcmpi(result.status,'OPTIMAL')
        v = result.x;
        sol.x_si = v(1:n);
        sol.z_si = v(n+1:n+m);
        sol.success = true;
    else
        sol.success = false;
        sol.x_si = [];
        sol.z_si = [];
    end
end
