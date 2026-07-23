function sol = solve_local_milp_intlinprog(A_si, b_si, EP_i, base0, x_si_k, g, params)
% Local MILP:
%   min   g' z
%   s.t.  A_si * x_si <= b_si
%         z >= base0 + EP_i*(x_si - x_si_k)
%         z >= 0
%         x_si in {0,1}^{nSi}
%
% Variables: v = [x_si; z]

    sol = struct('success', false, 'x_si', [], 'z', [], ...
                 'status', '', 'fval', [], 'exitflag', []);

    nSi = size(A_si,2);
    m   = length(base0);

    if length(g) ~= m
        error('Length of g must match length(base0) = #z.');
    end

    nvar = nSi + m;

    % objective
    f = [zeros(nSi,1); g(:)];

    % integer variables: x_si only
    intcon = 1:nSi;

    % bounds
    lb = [zeros(nSi,1); zeros(m,1)];
    ub = [ones(nSi,1);  inf(m,1)];

    % constraint (1): A_si x_si <= b_si
    A1 = sparse(size(A_si,1), nvar);
    A1(:,1:nSi) = sparse(A_si);
    b1 = b_si(:);

    % constraint (2):
    % z >= base0 + EP_i*(x_si - x_si_k)
    % -> EP_i*x_si - z <= -(base0 - EP_i*x_si_k)
    A2 = sparse(m, nvar);
    A2(:,1:nSi)     = sparse(EP_i);
    A2(:,nSi+1:end) = -speye(m);

    b2 = -(base0(:) - EP_i * x_si_k(:));

    A = [A1; A2];
    b = [b1; b2];

    % no equality constraints
    Aeq = [];
    beq = [];

    % options
    if nargin < 7 || isempty(params)
        options = optimoptions('intlinprog', ...
            'Display', 'off');
    else
        options = params;
    end

    x0 = [];

    try
       [v, fval, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, x0, options);
    catch ME
        sol.status = ['intlinprog_error: ', ME.message];
        return;
    end

    sol.exitflag = exitflag;
    sol.fval = fval;

    if exitflag > 0
        sol.success = true;
        sol.status = 'OPTIMAL_OR_FEASIBLE';
        sol.x_si = v(1:nSi);
        sol.z    = v(nSi+1:end);
    else
        sol.success = false;
        sol.status = output.message;
    end
end