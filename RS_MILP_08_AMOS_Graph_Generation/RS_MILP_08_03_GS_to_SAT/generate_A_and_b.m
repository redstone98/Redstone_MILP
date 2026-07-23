%% RS_MILP_07_GS_to_SAT
function [A, b, P_matrix, A_gj_info, b_gj_info, S_j_info, U_j_info, V_j_info] = generate_A_and_b(A_matrix, q, tau_vector, E1_Gj, E2_Gj_x, E2_Gj_t)
    A = [];
    b = [];
    P_matrix_transpose = [];
    
    
    N = length(A_matrix(:,1));

    U = 0;   
    for j = 1:q
      if isempty(E1_Gj.(['gs',num2str(j)])) 
      continue;
      end
    
      A_j = E2_Gj_x.(['gs',num2str(j)]) * E1_Gj.(['gs', num2str(j)]);
      A = [A;A_j];
      b_j_temp =  E2_Gj_t.(['gs',num2str(j)]) * E1_Gj.(['gs', num2str(j)]);
    
    
      A_gj_info.(['gs', num2str(j)]) = E2_Gj_x.(['gs', num2str(j)]);
      b_gj_info.(['gs', num2str(j)]) = ones(length(b_j_temp(:,1)),1) + 1/tau_vector(j) * b_j_temp* A_matrix(:,3);


      b = [b;b_gj_info.(['gs', num2str(j)])];
    
      P_matrix_transpose = [P_matrix_transpose; E1_Gj.(['gs', num2str(j)])];
    
      E1_gj = E1_Gj.(['gs', num2str(j)]);
      S_j_info.(['gs', num2str(j)]) = length(E1_gj(:,1));
    
      U_j_info.(['gs', num2str(j)]) = U;
      V_j_info.(['gs', num2str(j)]) = N - length(E1_gj(:,1)) - U;
      U = U + length(E1_gj(:,1));
    end
    
    
    P_matrix = P_matrix_transpose';

end