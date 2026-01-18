%% *RS-MILP-01: Configure Optimization Problem from Contact Matrix*
% 
% *University of Texas at Austin*
% 
% *Department of Aerospace Engineering and Engineering Mechanics*
% 
% *Hongseok Kim*
%% *1. Generate A matrix from contact chart*

clear;

addpath ~/Desktop/Redstone_MILP/RS_MILP_01_Config_MILP/
load('EOIR_48_SATs_12_Orbit_Planes_98_inc_7days_access_interval.mat','EOIR_access_interval')


% 기준 시각 정의

T = EOIR_access_interval;

t0 = datetime(2030, 1, 1, 0, 0, 0,'TimeZone','UTC');

% --- Source / Target을 string으로 통일 ---
src = string(T.Source);
tgt = string(T.Target);

% --- Ground 번호 추출 ---
gTok = regexp(src, 'Ground[_\s]*Point[_\s]*(\d+)', 'tokens', 'once');
groundNum = cellfun(@(c) str2double(c{1}), gTok);

% --- Satellite 번호 추출 ---
sTok = regexp(tgt, '(?i)SAT[^0-9]*(\d+)', 'tokens', 'once');
satNum = cellfun(@(c) str2double(c{1}), sTok);

% --- Start / End → datetime ---
st = T.StartTime;
en = T.EndTime;

if ~isdatetime(st)
    st = datetime(string(st), 'InputFormat','dd-MMM-yyyy HH:mm:ss');
end
if ~isdatetime(en)
    en = datetime(string(en), 'InputFormat','dd-MMM-yyyy HH:mm:ss');
end

% --- 중간 시각 ---
midTime = st + (en - st)/2;

% --- 기준 시각 대비 초(second) 차이 → integer ---
timeSec = seconds(midTime - t0);
timeSec = double(round(timeSec));   % 정수화 (필요시 floor / ceil로 변경 가능)




% --- 최종 3-column matrix ---
A_matrix = [double(satNum), double(groundNum), timeSec];
A_matrix = sortrows(A_matrix, 3);

A_matrix = A_matrix(1:100,:);
%
% 
% 
% 
%% 2. Extract Key parameters from A matrix

t = A_matrix(:,3);

% Satellite Cadence Constraint
tau = 60;

% Number of SATs
p = 48;

% Number of GSs
q = 54;

% Total number of contact
N = length(A_matrix(:,1));

% Number of contact for each SAT
S_i_vec = zeros(p,1);

for sat_index = 1:p
    S_i_vec(sat_index) = nnz(A_matrix(:,1) == sat_index);
end

% Number of contact for each GS
G_j_vec = zeros(q,1);
for gs_index = 1:q
    G_j_vec(gs_index) = nnz(A_matrix(:,2) == gs_index);
end
%% 
% 
%% 3. Selection matrix generation from given constant parameters
% 
% 3.1. $E_{S_i }^1$: Selection matrix from A-matrix to each satellite's contact sequence

E1_Si = struct();

for i = 1:p

E1_Si_mat = zeros(S_i_vec(i),N);

a_i = find(A_matrix(:,1) == i);

  for j = 1:length(a_i)
      E1_Si_mat(j, a_i(j)) = 1;
  end

  E1_Si.(['sat',num2str(i)]) = E1_Si_mat;

end
%% 
% 
% 
% 
% 3.2. $E_{S_i ,x}^2 \;,E_{S_{i,t} }^2$: Selection matrix for $\Delta t$ of each satelltite from $\left|S_i \right|$

E2_Si_x = struct();

for i = 1:p
    E2_Si_x_mat = [];
    for alpha = 1:S_i_vec(i)-1
        E2_Si_x_alpha = zeros(S_i_vec(i)-alpha, S_i_vec(i));
        E2_Si_x_alpha(:,alpha) = 1;
        for beta = alpha+1:S_i_vec(i)
            E2_Si_x_alpha(beta-alpha,beta) = 1;
        end
     E2_Si_x_mat = [E2_Si_x_mat;E2_Si_x_alpha]; 
    end
   if isempty(E2_Si_x_mat)
     E2_Si_x_mat = 0;
   end
   E2_Si_x.(['sat',num2str(i)]) = E2_Si_x_mat;
end


E2_Si_t = struct();

for i = 1:p
    E2_Si_t_mat = [];
    for alpha = 1:S_i_vec(i)-1
        E2_Si_t_alpha = zeros(S_i_vec(i)-alpha, S_i_vec(i));
        E2_Si_t_alpha(:,alpha) = -1;
        for beta = alpha+1:S_i_vec(i)
            E2_Si_t_alpha(beta-alpha,beta) = 1;
        end
     E2_Si_t_mat = [E2_Si_t_mat;E2_Si_t_alpha]; 
    end
   if isempty(E2_Si_t_mat)
     E2_Si_t_mat = 0;
   end
    E2_Si_t.(['sat',num2str(i)]) = E2_Si_t_mat;
end
%% 
% 
% 3.3. $E_{G_j }^1$: Selection matrix from A-matrix to each Ground Point's revisit sequence

E1_Gj = struct();

for j = 1:q
    E1_Gj_mat = zeros(G_j_vec(j)+2, N+2);

    E1_Gj_mat(1,1) = 1;
    b_j = find(A_matrix(:,2) == j);

    for alpha = 1:length(b_j)
        E1_Gj_mat(alpha+1, b_j(alpha)+1) = 1;
    end
    E1_Gj_mat(G_j_vec(j)+2, N+2) = 1;
    E1_Gj.(['gs', num2str(j)]) = E1_Gj_mat;
end
%% 
% 
% 3.4. $E_{G_{j,x} }^2 \;,E_{G_{j,t} }^2$: Selection matrix for $\Delta t$ of each GS from $\left|G_j \right|$

E2_Gj_x = struct();

for j = 1:q
    E2_Gj_x_mat = [];
    for alpha = 1:G_j_vec(j)+1
        E2_Gj_x_alpha = zeros(G_j_vec(j)+2-alpha, G_j_vec(j)+2);
        E2_Gj_x_alpha(:,alpha) = 1;
        for beta = alpha+1:G_j_vec(j)+2
            if alpha + 1 <= beta-1
               E2_Gj_x_alpha(beta-alpha, alpha+1:beta-1) = -1; 
            end
            E2_Gj_x_alpha(beta-alpha,beta) = 1;
        end
     E2_Gj_x_mat = [E2_Gj_x_mat;E2_Gj_x_alpha]; 
    end
   if isempty(E2_Gj_x_mat)
     E2_Gj_x_mat = 0;
   end
   E2_Gj_x.(['gs',num2str(j)]) = E2_Gj_x_mat;
end

E2_Gj_t = struct();

for j = 1:q
    E2_Gj_t_mat = [];
    for alpha = 1:G_j_vec(j)+1
        E2_Gj_t_alpha = zeros(G_j_vec(j)+2-alpha, G_j_vec(j)+2);
        E2_Gj_t_alpha(:,alpha) = -1;
        for beta = alpha+1:G_j_vec(j)+2
            E2_Gj_t_alpha(beta-alpha,beta) = 1;
        end
     E2_Gj_t_mat = [E2_Gj_t_mat;E2_Gj_t_alpha]; 
    end
   if isempty(E2_Gj_t_mat)
     E2_Gj_t_mat = 0;
   end
   E2_Gj_t.(['gs',num2str(j)]) = E2_Gj_t_mat;
end
%% 
% 
% 4. Derivation of $A,b,C,d,E,f,G$ revisit time problem to MILP
% 
% 4.1 $A,b$ for $L_1$ problem

A = [];
b_mat = [];
for i = 1:p
  if isempty(E1_Si.(['sat', num2str(i)])) 
  continue;
  end

  A_i = E2_Si_x.(['sat', num2str(i)]) * E1_Si.(['sat', num2str(i)]);
  A = [A;A_i];
  b_i_mat =  E2_Si_t.(['sat', num2str(i)]) * E1_Si.(['sat', num2str(i)]);
  b_mat = [b_mat; b_i_mat];
end

b = ones(length(b_mat(:,1)),1) + 1/tau * b_mat * A_matrix(:,3);

[A_row, A_col] = find(A==1);
figure;
scatter(A_row, A_col,'r','.')
grid on
figure;
scatter(1:length(b),b,'r','.')
grid on
%% 
% 
% 4.2 $C,d$ for $L_{\infty }$ problem

% C = [];
% d = [];
% 
% t = A_matrix(:,3);                  % (N×1) time vector in seconds (already prepared)
% t_aug = [0; t; max(A_matrix(:,3))+1];        % (N+2 × 1)
% 
% for j = 1:q
%     key = ['gs', num2str(j)];
% 
%     if isempty(E1_Gj.(key))
%         continue;
%     end
% 
%     % --- shorthand ---
%     E1  = E1_Gj.(key);        % selection matrix for this ground
%     E2x = E2_Gj_x.(key);
%     E2t = E2_Gj_t.(key);
% 
%     % =========================
%     % Build C_j  (constant)
%     % =========================
%     % D_j = diag( E2t*E1 * t_aug )   (vector -> diagonal matrix)
%     Dj_vec = (E2t * E1) * t_aug;           % (#rows_j × 1)
%     Dj = spdiags(Dj_vec, 0, length(Dj_vec), length(Dj_vec));  % sparse diag is safer
% 
%     % P = [0_{1×N}; I_{N×N}; 0_{1×N}]  -> (N+2 × N)
%     N = length(t);
%     P = [zeros(1,N); speye(N); zeros(1,N)];
% 
%     % C_j = D_j * (E2x*E1) * P
%     Cj = Dj * (E2x * E1) * P;
% 
%     % =========================
%     % Build d_j  (constant)
%     % =========================
%     % v = [1; 0_{N×1}; 1] - 1  = [0; -1...; 0]
%     v = [1; zeros(N,1); 1];          % (N+2 × 1)
% 
%     % d_j = D_j * (E2x*E1) * v
%     dj = Dj * ((E2x * E1) * v - ones(length(Dj(:,1)),1));
% 
%     % --- stack ---
%     C = [C; Cj];
%     d = [d; dj];
% end
% 
% 
% [C_row, C_col] = find(C~=0);
% 
% figure;
% scatter(C_row, C_col,'r','.')
% grid on
%% 
% 
% 
% 
% 4.3 $E,f,G$ for $L_2$ Optimization problem

% % =========================
% % Build E, f, G  (2nd image)
% % =========================
% 
% E = [];
% F = [];
% 
% t = A_matrix(:,3);                       % (N×1)
% N = length(t);
% 
% t_aug = [0; t; max(t)+1];                % (N+2 × 1), same as your code
% 
% % P = [0_{1×N}; I_N; 0_{1×N}]  -> (N+2 × N)
% P = [zeros(1,N); speye(N); zeros(1,N)];
% 
% % this is the (N+2 × 1) vector [1; 0_N; 1]
% one_aug = [1; zeros(N,1); 1];
% 
% % For G, we will build the stacked matrix Bt = [E2t*E1; ...] first (optional),
% % or directly build gvec = Bt*t_aug.
% gvec = [];   % = Bt * t_aug  (will be stacked across j)
% 
% for j = 1:q
%     key = ['gs', num2str(j)];
% 
%     if isempty(E1_Gj.(key))
%         continue;
%     end
% 
%     % --- shorthand ---
%     E1  = E1_Gj.(key);        % selection matrix for this ground
%     E2x = E2_Gj_x.(key);
%     E2t = E2_Gj_t.(key);
% 
%     % Common blocks
%     Sx = (E2x * E1);          % (#rows_j × (N+2))
%     St = (E2t * E1);          % (#rows_j × (N+2))
% 
%     % -------------------------
%     % E_j = Sx * P    (constant)
%     % -------------------------
%     Ej = Sx;              % (#rows_j × N)
%     E  = [E; Ej];
% 
%     % -------------------------
%     % f_j = Sx]
%     % -------------------------
%     Fj = Sx;        % (#rows_j × 1)
%     F  = [Fj; Fj];
% 
%     % -------------------------
%     % gvec_j = St * t_aug
%     % -------------------------
%     gvec_j = St * t_aug;      % (#rows_j × 1)
%     gvec   = [gvec; gvec_j];
% end
% 
% 
% E = E * P;
% f = F * one_aug;
% 
% % -------------------------
% % G = gvec * gvec'
% % -------------------------
% % WARNING: this can be very large/dense if gvec is long.
% % If you only need y'*G*y, note that y'*(gvec*gvec')*y = (gvec'*y)^2.
% G = gvec * gvec.';   % (M×M)
% 
% 
% [E_row, E_col] = find(E~=0);
% 
% figure;
% scatter(E_row, E_col,'r','.')
% grid on
% imagesc(G); colorbar; axis equal tight;
% title('Heatmap of G');
%% 
% 
% 
% 
%% 5. Optimization Problem (MILP)

L1_flag = 1;
L_infty_flag = 0;
%% 
% 
% 5.1. $L_1$ revisit time problem to MILP

if L1_flag == 1

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

[x_opt, fval, exitflag, output] = intlinprog( ...
    f, intcon, Aineq, bineq, [], [], lb, ub);
% ---------- 5) Recover maximization objective ----------
max_onesTx = -fval;       % because f = -1
x_opt = round(x_opt);     % safety: should already be integer

% ---------- 6) Quick checks ----------
viol = Aineq*x_opt - bineq;
max_viol = max(viol);

fprintf('Exitflag: %d\n', exitflag);
fprintf('Objective (max 1^T x): %.0f\n', max_onesTx);
fprintf('Max constraint violation: %.3e\n', max_viol);
end
%% 
% 
% 5.2 $L_{\infty }$ revisit time problem to MILP

if L_infty_flag == 1

%% Given:
% A (m×N), b (m×1)
% C (q×N), d (q×1)
% N (scalar)

% ---------- sanity ----------
[m, nA] = size(A);
assert(nA == N, 'A must have N columns.');
b = b(:);  assert(length(b) == m, 'b must be m×1.');

[q, nC] = size(C);
assert(nC == N, 'C must have N columns.');
d = d(:);  assert(length(d) == q, 'd must be q×1.');

% ---------- decision variables ----------
% z = [x; R]  -> length N+1
nvar = N + 1;

% Objective: min R  => f = [0...0, 1]
f = [zeros(N,1); 1];

% Integer constraints: x binary, R continuous
intcon = 1:N;

% Bounds
lb = [zeros(N,1); 0];      % R >= 0 (필요 없으면 -Inf로 바꿔도 됨)
ub = [ones(N,1);  Inf];

% ---------- Build inequalities Aineq*z <= bineq ----------
% 1) Ax <= b  -> [A, 0] [x;R] <= b
A1 = [A, zeros(m,1)];
b1 = b;

% 2) Cx - R*1 <= -d  -> [C, -1] [x;R] <= -d   (각 row마다 -R)
A2 = [C, -ones(q,1)];
b2 = -d;

Aineq = [A1; A2];
bineq = [b1; b2];

% (No equality constraints)
Aeq = [];
beq = [];

% ---------- Solve ----------
% opts = optimoptions('intlinprog', ...
%     'Display','iter', ...
%     'Heuristics','advanced', ...
%     'CutGeneration','advanced');

[z_opt, fval, exitflag, output] = intlinprog( ...
    f, intcon, Aineq, bineq, Aeq, beq, lb, ub);
% ---------- Parse solution ----------
x_opt = round(z_opt(1:N));     % binary
R_opt = z_opt(N+1);            % optimal R

% ---------- Feasibility check ----------
viol1 = A*x_opt - b;
viol2 = C*x_opt + d - R_opt*ones(q,1);   % should be <= 0
fprintf('Exitflag: %d\n', exitflag);
fprintf('R_opt: %.6f\n', R_opt);
fprintf('max(Ax-b): %.3e\n', max(viol1));
fprintf('max(Cx+d-R): %.3e\n', max(viol2));
end
%% 
% 
% 
% 
% 
%