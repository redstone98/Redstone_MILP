function sol = solve_local_milp_gurobi(A_si, b_si, EP_i, base0, x_si_k, g, params)
% Local MILP:
%   min   g' z
%   s.t.  A_si * x_si <= b_si
%         z >= base0 + EP_i*(x_si - x_si_k)
%         z >= 0
%         x_si in {0,1}^{nSi}
%
% Variables: v = [x_si; z]

    sol = struct('success', false, 'x_si', [], 'z', [], 'status', '');

    nSi = size(A_si,2);
    m   = length(base0);

    if length(g) ~= m
        error('Length of g must match length(base0) = #z.');
    end

    nvar = nSi + m;
    model.modelsense = 'min';

    % ----- linear objective: g' z -----
    model.obj = [zeros(nSi,1); g(:)];

    % ----- variable types -----
    model.vtype = [repmat('B', nSi, 1); repmat('C', m, 1)];

    % ----- bounds -----
    model.lb = [zeros(nSi,1); zeros(m,1)];
    model.ub = [ones(nSi,1);  inf(m,1)];

    % ----------------------------
    % Constraints
    % (1) A_si x_si <= b_si
    % (2) z >= base0 + EP_i*(x_si - x_si_k)
    %     -> -z + EP_i*x_si <= -(base0 - EP_i*x_si_k)
    % ----------------------------
    % (1)
    A1 = sparse(size(A_si,1), nvar);
    A1(:,1:nSi) = sparse(A_si);
    rhs1 = b_si(:);
    sense1 = repmat('<', size(A_si,1), 1);

    % (2)
    A2 = sparse(m, nvar);
    A2(:,1:nSi)     = sparse(EP_i);
    A2(:,nSi+1:end) = -speye(m);

    rhs2 = -(base0(:) - EP_i * x_si_k(:));
    sense2 = repmat('<', m, 1);

    model.A     = [A1; A2];
    model.rhs   = [rhs1; rhs2];
    model.sense = [sense1; sense2];

    % params
    if nargin < 7 || isempty(params)
        params = struct();
    end

    % solve
    try
        result = gurobi(model, params);
    catch ME
        sol.status = ['gurobi_error: ', ME.message];
        return;
    end

    sol.status = result.status;

    if strcmpi(result.status,'OPTIMAL') || strcmpi(result.status,'SUBOPTIMAL')
        v = result.x;
        sol.x_si = v(1:nSi);
        sol.z    = v(nSi+1:end);
        sol.success = true;
    end
end