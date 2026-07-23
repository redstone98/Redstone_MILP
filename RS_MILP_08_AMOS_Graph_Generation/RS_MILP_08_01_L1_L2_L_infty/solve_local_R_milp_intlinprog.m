function sol = solve_local_R_milp_intlinprog( ...
    A_si, b_si, EP_i, base0, x_si_k, gvec, options)
% Local BCD MILP using MATLAB intlinprog
%
%   min_{x_si, R_si} R_si
%
%   s.t.
%       A_si*x_si <= b_si
%
%       R_si*1 >= diag(gvec) * ...
%           (base0 + EP_i*(x_si - x_si_k))
%
% where:
%   x_si in {0,1}^nSi
%   R_si >= 0
%
% Variables:
%   v = [x_si; R_si]

    sol = struct( ...
        'success', false, ...
        'x_si', [], ...
        'R', [], ...
        'status', '', ...
        'exitflag', [], ...
        'output', [], ...
        'objective', []);

    if nargin < 7 || isempty(options)
        options = optimoptions('intlinprog', ...
            'Display', 'off');
    end

    % =====================================================
    % Preprocessing
    % =====================================================
    A_si   = sparse(A_si);
    EP_i   = sparse(EP_i);

    b_si   = b_si(:);
    base0  = base0(:);
    x_si_k = x_si_k(:);
    gvec   = gvec(:);

    nSi = size(A_si, 2);
    m   = length(base0);

    if size(EP_i, 1) ~= m
        error('size(EP_i,1) must match length(base0).');
    end

    if size(EP_i, 2) ~= nSi
        error('size(EP_i,2) must match size(A_si,2).');
    end

    if length(x_si_k) ~= nSi
        error('length(x_si_k) must match size(A_si,2).');
    end

    if length(gvec) ~= m
        error('length(gvec) must match length(base0).');
    end

    if length(b_si) ~= size(A_si,1)
        error('length(b_si) must match size(A_si,1).');
    end

    % =====================================================
    % Variables:
    %
    % v = [x_si; R_si]
    % =====================================================
    nvar = nSi + 1;

    % Objective:
    % min R_si
    f = [zeros(nSi,1); 1];

    % Binary/integer variable indices
    intcon = 1:nSi;

    % Bounds
    lb = [zeros(nSi,1); 0];
    ub = [ones(nSi,1); inf];

    % =====================================================
    % Constraint 1:
    %
    % A_si*x_si <= b_si
    % =====================================================
    A1 = [A_si, sparse(size(A_si,1),1)];
    rhs1 = b_si;

    % =====================================================
    % Constraint 2:
    %
    % R_si*1 >= diag(gvec) * ...
    %     (base0 + EP_i*(x_si - x_si_k))
    %
    % Define:
    %
    % base_i = base0 - EP_i*x_si_k
    %
    % Then:
    %
    % R_si*1 >= diag(gvec)*(base_i + EP_i*x_si)
    %
    % Equivalent <= form:
    %
    % diag(gvec)*EP_i*x_si - R_si*1
    %     <= -diag(gvec)*base_i
    % =====================================================
    base_i = base0 - EP_i*x_si_k;

    G = spdiags(gvec, 0, m, m);

    A2 = [G*EP_i, -sparse(ones(m,1))];
    rhs2 = -gvec .* base_i;

    % Combined inequality constraints
    Aineq = sparse([A1; A2]);
    bineq = [rhs1; rhs2];

    % No equality constraints
    Aeq = [];
    beq = [];

    % =====================================================
    % Solve using intlinprog
    % =====================================================
    try
        [v, fval, exitflag, output] = intlinprog( ...
            f, intcon, ...
            Aineq, bineq, ...
            Aeq, beq, ...
            lb, ub, ...
            options);

    catch ME
        sol.status = ['intlinprog_error: ', ME.message];
        return;
    end

    % =====================================================
    % Store results
    % =====================================================
    sol.exitflag = exitflag;
    sol.output = output;

    if isempty(v)
        sol.status = get_intlinprog_status(exitflag);
        return;
    end

    sol.x_si = full(double(v(1:nSi)));
    sol.R = full(double(v(end)));
    sol.objective = full(double(fval));
    sol.status = get_intlinprog_status(exitflag);

    % exitflag > 0 indicates a successful solution
    sol.success = exitflag > 0;
end


function status = get_intlinprog_status(exitflag)
% Convert intlinprog exit flag to readable status text.

    switch exitflag
        case 1
            status = 'OPTIMAL';

        case 2
            status = 'FEASIBLE_INTEGER_SOLUTION';

        case 0
            status = 'LIMIT_REACHED';

        case -1
            status = 'STOPPED_BY_OUTPUT_FUNCTION';

        case -2
            status = 'INFEASIBLE';

        case -3
            status = 'UNBOUNDED_OR_ROOT_LP_UNBOUNDED';

        case -9
            status = 'SOLVER_ERROR';

        otherwise
            status = sprintf('UNKNOWN_EXITFLAG_%d', exitflag);
    end
end