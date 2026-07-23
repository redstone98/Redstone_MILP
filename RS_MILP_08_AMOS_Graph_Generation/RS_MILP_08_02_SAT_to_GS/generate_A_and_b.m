function [A, b, P_matrix, A_si_info, b_si_info, S_i_info, U_i_info, V_i_info] = generate_A_and_b(A_matrix, p, tau, E1_Si, E2_Si_x, E2_Si_t)
    A = [];
    b_mat = [];
    P_matrix_transpose = [];
    
    
    N = length(A_matrix(:,1));

    U = 0;   
    for i = 1:p
      if isempty(E1_Si.(['sat', num2str(i)])) 
      continue;
      end
    
      A_i = E2_Si_x.(['sat', num2str(i)]) * E1_Si.(['sat', num2str(i)]);
      A = [A;A_i];
      b_i_mat =  E2_Si_t.(['sat', num2str(i)]) * E1_Si.(['sat', num2str(i)]);
      b_mat = [b_mat; b_i_mat];
    
    
      A_si_info.(['sat', num2str(i)]) = E2_Si_x.(['sat', num2str(i)]);
      b_si_info.(['sat', num2str(i)]) = ones(length(b_i_mat(:,1)),1) + 1/tau * b_i_mat* A_matrix(:,3);
    
      P_matrix_transpose = [P_matrix_transpose; E1_Si.(['sat', num2str(i)])];
    
      E1_si = E1_Si.(['sat', num2str(i)]);
      S_i_info.(['sat', num2str(i)]) = length(E1_si(:,1));
    
      U_i_info.(['sat', num2str(i)]) = U;
      V_i_info.(['sat', num2str(i)]) = N - length(E1_si(:,1)) - U;
      U = U + length(E1_si(:,1));
    
    end
    
    
    
    b = ones(length(b_mat(:,1)),1) + 1/tau * b_mat * A_matrix(:,3);
    P_matrix = P_matrix_transpose';

end