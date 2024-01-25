function K = K_fn(m, n)

    % Get commutation matrix (Lutkepohl (2005, Appendix A12.2))
    
    % Inputs:
    % m         1 x 1                 commutation dimension #1
    % n         1 x 1                 commutation dimension #2
    
    % Output:
    % D_pl     (m*n) x (m*n)          K matrix
    

    K = zeros(m*n);
    I = eye(m*n);
    k = 1;
    for i=1:n
        for j=i:n:m*n
            K(k,:) = I(j,:);
            k = k+1;
        end
    end

end