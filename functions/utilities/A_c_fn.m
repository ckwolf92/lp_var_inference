function A_c = A_c_fn(A)

    % Get companion form matrix
    
    % Inputs:
    % A         n_y x (n_y*p)     matrix of VAR coefficients (without constant)
    
    % Output:
    % A_c       (n_y*p) x (n_y*p) companion for matrix
    

    n = size(A,1);
    p = size(A,2)/n;
    
    if p == 1
        A_c = A;
    else
        A_c = zeros(p*n,p*n);
        A_c(1:n,:) = A;
        A_c(n+1:n*p,1:n*(p-1)) = eye(n*(p-1));
    end

end