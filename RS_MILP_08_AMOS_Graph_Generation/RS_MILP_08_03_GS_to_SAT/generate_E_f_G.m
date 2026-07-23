function [E, f_vector, gvec] = generate_E_f_G(A_matrix, p, E1_Si, E2_Si_x, E2_Si_t, t_start, t_end, k_vector)

E_ori = [];
F = [];

t = A_matrix(:,3);                       % (N×1)
N = length(t);

t_aug = [t_start; t; t_end];                % (N+2 × 1), same as your code

% P = [0_{1×N}; I_N; 0_{1×N}]  -> (N+2 × N)
P = [zeros(1,N); speye(N); zeros(1,N)];

% this is the (N+2 × 1) vector [1; 0_N; 1]
one_aug = [1; zeros(N,1); 1];

% For G, we will build the stacked matrix Bt = [E2t*E1; ...] first (optional),
% or directly build gvec = Bt*t_aug.
gvec = [];   % = Bt * t_aug  (will be stacked across j)

for i = 1:p
    key = ['sat', num2str(i)];

    if isempty(E1_Si.(key))
        continue;
    end

    % --- shorthand ---
    E1  = E1_Si.(key);        % selection matrix for this ground
    E2x = E2_Si_x.(key);
    E2t = E2_Si_t.(key);

    % Common blocks
    Sx = (E2x * E1);          % (#rows_j × (N+2))
    St = (E2t * E1);          % (#rows_j × (N+2))

    % -------------------------
    % E_j = Sx * P    (constant)
    % -------------------------
    Ej = Sx;              % (#rows_j × N)
    E_ori  = [E_ori; Ej];

    % -------------------------
    % f_j = Sx]
    % -------------------------
    Fj = Sx;        % (#rows_j × 1)
    F  = [F; Fj];

    % -------------------------
    % gvec_j = St * t_aug
    % -------------------------
    gvec_j = k_vector(i) * St * t_aug;      % (#rows_j × 1)
    gvec   = [gvec; gvec_j];
end

% Dj = spdiasat(Dj_vec, 0, length(Dj_vec), length(Dj_vec)); 

E = E_ori * P;
f_vector = F * one_aug;

end