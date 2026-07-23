function  x_L_infinity = solve_L_infty(x, N, p, A_si_info, b_si_info, E, P_matrix, f_vector, gvec ,E1_Si, S_i_info, U_i_info)
    
    % ---- user provided initializations (keep) ----
    x_history_BCD = [];
    R_history_BCD = [];
    
    % ---- BCD options ----
    maxOuterIters = 50;          % outer sweeps (full i=1..p)
    tol_x_change  = 0;           % since binary, exact equality check is fine
    verbose       = 1;
    
    % ---- precompute constant matrices ----
    % P_matrix: permutation such that x = P * [x_s1; ...; x_sp]
    % EP = E * P  (used in all agents)
    EP = E * P_matrix;
    
    % ---- build block index map (idx_i for each agent) ----
    % You already have:
    %   U = sum_{j<i} |S_j|,  S_i is list (or just length defines block size)
    % We'll build idx_i = (U+1):(U+|S_i|)
    idx_block = cell(p,1);
    for i = 1:p
        if isempty(E1_Si.(['sat', num2str(i)])) 
            idx_block{i} = [];
            continue;
        end
        U = U_i_info.(['sat', num2str(i)]);
        S_i = S_i_info.(['sat', num2str(i)]);
        nSi = length(S_i);
        idx_block{i} = (U+1):(U+nSi);
    end
    
    BCD_handle = 1;
    outerIter  = 0;
    
    while BCD_handle == 1
        outerIter = outerIter + 1;
        if outerIter > maxOuterIters
            if verbose, fprintf('[BCD] Reached maxOuterIters=%d\n', maxOuterIters); end
            break;
        end
    
        x_old_sweep = x;
    
        for i = 1:p
            % skip empty agent
            if isempty(E1_Si.(['sat', num2str(i)]))
                continue;
            end
    
            % ----- snapshot current iterate -----
            x_k = x;
    
            % ----- local constraints -----
            A_si = A_si_info.(['sat', num2str(i)]);
            b_si = b_si_info.(['sat', num2str(i)]);
    
            % ----- stacked vector at k: xbar^k = P' * x^k -----
            xbar_k = P_matrix' * x_k;
    
            % ----- identify current local block -----
            idx_i  = idx_block{i};
            x_si_k = E1_Si.(['sat', num2str(i)]) * x_k;       % this is the "current" local block
    
            % ----- coupling base term: base0 = EP*xbar_k + f - 1 -----
            base0  = EP * xbar_k + f_vector - 1;
    
            % ----- the only part that multiplies x_si variable -----
    
    
            U = U_i_info.(['sat', num2str(i)]);
            S_i = S_i_info.(['sat', num2str(i)]);
            EP_i = EP(:, (U+1):(U+S_i));          % == EP * J_i
    
            J_i = sparse(N, S_i);              % (U+|S_i|+V) x |S_i|
            J_i(U+1:U+S_i, :) = speye(S_i);    % middle block = I_|S_i|
    
            % solve local subproblem with Gurobi MIQP
            params = struct();
            params.OutputFlag = 0;   % set 1 for debugging
    %%%%%
    
            sol = solve_local_R_milp_gurobi(A_si, b_si, EP_i, base0, x_si_k, gvec, params);
            
            if ~sol.success
                fprintf('[BCD] sat%d local solve failed: %s\n', i, sol.status);
                continue;
            end
            
            x_si_new = sol.x_si;
            R_new = sol.R;
            
            xbar_new = xbar_k + J_i*(x_si_new - x_si_k);
            x = P_matrix * xbar_new;
            
            R_history_BCD = [R_history_BCD, R_new];
            x_history_BCD = [x_history_BCD, x];
            x_L_infinity = x;
     %%%%
        end
    
        % ----- stopping check (end of sweep) -----
        if isequal(x, x_old_sweep) <= 1e-12
            if verbose
                fprintf('[BCD] Converged at outerIter=%d (no change after full sweep)\n', outerIter);
            end
            break;
        else
            if verbose
                dx = nnz(x ~= x_old_sweep);
                fprintf('[BCD] outerIter=%d done. nnz(dx)=%d', outerIter, dx);
            end
        end
    end

end