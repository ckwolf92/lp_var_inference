function H_c = H_c_fn(H,p)

    % Get companion form H matrix
    
    % Inputs:
    % H         n_y x n_y       H matrix of local-to-VAR model (structural shock impact IRF)
    
    % Output:
    % H_c       (n_y*p) x n_y   companion form of H
    

    n   = size(H,1);
    H_c = [H;zeros((p-1)*n,n)];

end