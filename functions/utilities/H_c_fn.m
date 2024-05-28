function H_c = H_c_fn(H,p)

    % Get companion form matrix
    
    % Inputs:
    % A         n_y x (n_y*p)     matrix of VAR coefficients (without constant)
    
    % Output:
    % A_c       (n_y*p) x (n_y*p) companion form matrix
    

    n   = size(H,1);
    H_c = [H;zeros((p-1)*n,n)];

end