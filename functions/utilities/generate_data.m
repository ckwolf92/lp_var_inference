function data_y = generate_data(dgp)

    % Simulate data ABCD system
    
    % Inputs:
    % A         n_s x n_s       state evolution
    % B         n_s x n_e       state shock response
    % A         n_y x n_s       observables evolution
    % D         n_y x n_e       observables shock response
    % T         1 x 1           sample size
    
    % Output:
    % data_y    (T+1) x 1       simulated data (with initial condition)


    % preparations
    
    T = dgp.T;
    
    A = dgp.ABCD.A;
    B = dgp.ABCD.B;
    C = dgp.ABCD.C;
    D = dgp.ABCD.D;
    
    n_s   = dgp.n_s;
    n_eps = dgp.n_eps;
    n_y   = dgp.n_y;
    
    % draw shocks
    
    data_eps = randn(T,n_eps);
    
    % simulate states
    
    s = zeros(n_s,1);
    data_s = NaN(T,n_s);
    for t = 1:T
        s = A * s + B * data_eps(t,:)';
        data_s(t,:) = s';
    end
    
    % simulate observables

    data_y = NaN(T,n_y);
    data_y(1,:) = (D * data_eps(1,:)')';
    for t = 2:T
        data_y(t,:) = (C * data_s(t-1,:)' + D * data_eps(t,:)')';
    end

    data_y = [zeros(1,n_y);data_y];

end