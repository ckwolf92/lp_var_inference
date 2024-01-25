function L = L_fn(n)

    % Get elimination matrix (Lutkepohl (2005, Appendix A12.2))
    
    % Inputs:
    % n         1 x 1                 system size
    
    % Output:
    % L         (n*(n+1)/2) x n^2     L matrix
    

    T = tril(ones(n));
    f = find(T(:));
    k = n*(n+1)/2;
    n2 = n*n;
    L = zeros(n2,k);
    x = f + n2*(0:k-1)';
    L(x) = 1;
    L = L';

end