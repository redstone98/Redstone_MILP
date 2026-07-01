function x_L1 = solve_L1(A,b,N)
 %% Given: A (m×N), b (m×1), N (scalar)
% Goal: max 1^T x  s.t. A x <= b, x in {0,1}^N
% ---------- 1) Dimensions / sanity checks ----------
[m, nA] = size(A);
assert(nA == N, 'A must have N columns.');
assert(isvector(b) && length(b) == m, 'b must be m×1 to match A.');

b = b(:);                 % force column vector
Aineq = A;
bineq = b;

% ---------- 2) Convert max to min ----------
% intlinprog solves: min f'*x
f = -ones(N,1);           % minimize -sum(x) == maximize sum(x)

% ---------- 3) Binary variable settings ----------
intcon = 1:N;             % all variables are integer
lb = zeros(N,1);
ub = ones(N,1);

% (Optional) If you truly want "binary", keep intcon + bounds [0,1].
% Alternatively, you can also set intcon and bounds; that's standard.

% ---------- 4) Solve MILP ----------
% opts = optimoptions('intlinprog', ...
%     'Display','iter', ...
%     'Heuristics','advanced', ...
%     'CutGeneration','advanced');

[x_opt_L1, fval, exitflag, output] = intlinprog( ...
    f, intcon, Aineq, bineq, [], [], lb, ub);
% ---------- 5) Recover maximization objective ----------
max_onesTx = -fval;       % because f = -1
x_L1 = round(x_opt_L1);     % safety: should already be integer

% ---------- 6) Quick checks ----------
viol = Aineq*x_opt_L1 - bineq;
max_viol = max(viol);

fprintf('Exitflag: %d\n', exitflag);
fprintf('Objective (max 1^T x): %.0f\n', max_onesTx);
fprintf('Max constraint violation: %.3e\n', max_viol);
end