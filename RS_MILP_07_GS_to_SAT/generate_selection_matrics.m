%% RS_MILP_07_GS_to_SAT


function [E1_Si, E2_Si_x, E2_Si_t, E1_Gj, E2_Gj_x, E2_Gj_t] = generate_selection_matrics(A_matrix, p,q)

%% 3.0 Get number of contacts for each SAT/GS

N = size(A_matrix,1);

% Number of contacts for each SAT
S_i_vec = zeros(p,1);
for i = 1:p
    S_i_vec(i) = nnz(A_matrix(:,1) == i);
end

% Number of contacts for each GS
G_j_vec = zeros(q,1);
for j = 1:q
    G_j_vec(j) = nnz(A_matrix(:,2) == j);
end


%% 3.1 E1_Gj: Selection matrix from A-matrix to each GS contact sequence

E1_Gj = struct();

for j = 1:q
    E1_Gj_mat = zeros(G_j_vec(j), N);
    b_j = find(A_matrix(:,2) == j);
    for alpha = 1:length(b_j)
        E1_Gj_mat(alpha, b_j(alpha)) = 1;
    end
    E1_Gj.(['gs', num2str(j)]) = E1_Gj_mat;

end


%% 3.2 E2_Gj_x, E2_Gj_t: Selection matrix for delta_t of each GS

E2_Gj_x = struct();

for j = 1:q
    E2_Gj_x_mat = [];
    for alpha = 1:G_j_vec(j)-1
        E2_Gj_x_alpha = zeros(G_j_vec(j)-alpha, G_j_vec(j));
        E2_Gj_x_alpha(:,alpha) = 1;
        for beta = alpha+1:G_j_vec(j)
            E2_Gj_x_alpha(beta-alpha,beta) = 1;
        end
        E2_Gj_x_mat = [E2_Gj_x_mat; E2_Gj_x_alpha];
    end

    if isempty(E2_Gj_x_mat)
        E2_Gj_x_mat = 0;
    end

    E2_Gj_x.(['gs', num2str(j)]) = E2_Gj_x_mat;

end


E2_Gj_t = struct();

for j = 1:q
    E2_Gj_t_mat = [];
    for alpha = 1:G_j_vec(j)-1
        E2_Gj_t_alpha = zeros(G_j_vec(j)-alpha, G_j_vec(j));
        E2_Gj_t_alpha(:,alpha) = -1;
        for beta = alpha+1:G_j_vec(j)
            E2_Gj_t_alpha(beta-alpha,beta) = 1;
        end
        E2_Gj_t_mat = [E2_Gj_t_mat; E2_Gj_t_alpha];
    end

    if isempty(E2_Gj_t_mat)
        E2_Gj_t_mat = 0;
    end

    E2_Gj_t.(['gs', num2str(j)]) = E2_Gj_t_mat;

end


%% 3.3 E1_Si: Selection matrix from A-matrix to each SAT contact sequence with boundaries

E1_Si = struct();

for i = 1:p
    E1_Si_mat = zeros(S_i_vec(i)+2, N+2);
    % start boundary
    E1_Si_mat(1,1) = 1;
    a_i = find(A_matrix(:,1) == i);
    for alpha = 1:length(a_i)
        E1_Si_mat(alpha+1, a_i(alpha)+1) = 1;
    end
    % end boundary
    E1_Si_mat(S_i_vec(i)+2, N+2) = 1;
    E1_Si.(['sat', num2str(i)]) = E1_Si_mat;
end


%% 3.4 E2_Si_x, E2_Si_t: Matrix for delta_t of each SAT with boundaries

E2_Si_x = struct();

for i = 1:p
    E2_Si_x_mat = [];
    for alpha = 1:S_i_vec(i)+1
        E2_Si_x_alpha = zeros(S_i_vec(i)+2-alpha, S_i_vec(i)+2);
        E2_Si_x_alpha(:,alpha) = 1;
        for beta = alpha+1:S_i_vec(i)+2
            if alpha + 1 <= beta - 1
                E2_Si_x_alpha(beta-alpha, alpha+1:beta-1) = -1;
            end
            E2_Si_x_alpha(beta-alpha,beta) = 1;
        end
        E2_Si_x_mat = [E2_Si_x_mat; E2_Si_x_alpha];
    end

    if isempty(E2_Si_x_mat)
        E2_Si_x_mat = 0;
    end

    E2_Si_x.(['sat', num2str(i)]) = E2_Si_x_mat;

end


E2_Si_t = struct();

for i = 1:p
    E2_Si_t_mat = [];
    for alpha = 1:S_i_vec(i)+1
        E2_Si_t_alpha = zeros(S_i_vec(i)+2-alpha, S_i_vec(i)+2);
        E2_Si_t_alpha(:,alpha) = -1;

        for beta = alpha+1:S_i_vec(i)+2
            E2_Si_t_alpha(beta-alpha,beta) = 1;
        end
        E2_Si_t_mat = [E2_Si_t_mat; E2_Si_t_alpha];
    end

    if isempty(E2_Si_t_mat)
        E2_Si_t_mat = 0;
    end

    E2_Si_t.(['sat', num2str(i)]) = E2_Si_t_mat;
end

end