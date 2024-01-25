function D_pl = D_pl_fn(n)

    % Get D+ matrix (Lutkepohl (2005, Appendix A12.2))
    
    % Inputs:
    % n         1 x 1                 system size
    
    % Output:
    % D_pl     (n*(n+1)/2) x n^2      D+ matrix
    

    m   = n * (n + 1) / 2;
    nsq = n^2;
    r   = 1;
    a   = 1;
    v   = zeros(1, nsq);
    cn  = cumsum(n:-1:2);
    for i = 1:n
       v(r:r + i - 2) = i - n + cn(1:i - 1);
       r = r + i - 1;
       
       v(r:r + n - i) = a:a + n - i;
       r = r + n - i + 1;
       a = a + n - i + 1;
    end
    D2 = sparse(1:nsq, v, 1, nsq, m);
    
    D_pl = pinv(full(D2));

end