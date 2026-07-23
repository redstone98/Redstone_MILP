%% RS_MILP_07_GS_to_SAT
function sol = solve_local_milp_intlinprog(A_gj, b_gj, EP_j, base0, x_gj_k, g, params)
% Local Gj-based MILP:
%
%   min   g' z
%
%   s.t.  A_gj * x_gj <= b_gj
%
%         z >= base0 + EP_j*(x_gj - x_gj_k)
%
%         z >= 0
%
%         x_gj in {0,1}^{nGj}
%
% Variables:
%   v = [x_gj; z]

    sol = struct('success', false, 'x_gj', [], 'z', [], ...
                 'status', '', 'fval', [], 'exitflag', []);

    nGj = size(A_gj,2);
    m   = length(base0);

    if length(g) ~= m
        error('Length of g must match length(base0) = #z.');
    end

    nvar = nGj + m;

    % objective:
    %
    %   min g' z
    %
    % Since v = [x_gj; z], the objective vector is:
    %
    %   f = [0; g]
    %
    f = [zeros(nGj,1); g(:)];

    % integer variables:
    % only x_gj is binary/integer
    intcon = 1:nGj;

    % bounds:
    %
    %   0 <= x_gj <= 1
    %   z >= 0
    %
    lb = [zeros(nGj,1); zeros(m,1)];
    ub = [ones(nGj,1);  inf(m,1)];

    % constraint (1):
    %
    %   A_gj x_gj <= b_gj
    %
    A1 = sparse(size(A_gj,1), nvar);
    A1(:,1:nGj) = sparse(A_gj);
    b1 = b_gj(:);

    % constraint (2):
    %
    %   z >= base0 + EP_j*(x_gj - x_gj_k)
    %
    % equivalently:
    %
    %   z >= base0 + EP_j*x_gj - EP_j*x_gj_k
    %
    %   EP_j*x_gj - z <= -base0 + EP_j*x_gj_k
    %
    %   EP_j*x_gj - z <= -(base0 - EP_j*x_gj_k)
    %
    A2 = sparse(m, nvar);
    A2(:,1:nGj)     = sparse(EP_j);
    A2(:,nGj+1:end) = -speye(m);

    b2 = -(base0(:) - EP_j * x_gj_k(:));

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
       [v, fval, exitflag, output] = intlinprog( ...
           f, intcon, A, b, Aeq, beq, lb, ub, x0, options);
    catch ME
        sol.status = ['intlinprog_error: ', ME.message];
        return;
    end

    sol.exitflag = exitflag;
    sol.fval = fval;

    if exitflag > 0
        sol.success = true;
        sol.status = 'OPTIMAL_OR_FEASIBLE';

        sol.x_gj = v(1:nGj);
        sol.z    = v(nGj+1:end);
    else
        sol.success = false;
        sol.status = output.message;
    end
end