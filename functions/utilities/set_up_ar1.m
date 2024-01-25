% Bring AR model into state-space form & compute IRFs

    % Original system, AR(1)
    % y_t = \rho Y_{t-1} + (1 + T^{-\zeta} * \alpha(L)) * \epsilon_t
    
    % State-space form
    % s_t = A s_{t-1} + B \epsilon_t
    % x_t = C s_{t-1} + D \epsilon_t
    % where
    % s_t       (y_t, \epsilon_t, \epsilon_{t-1}, ...)'
    % x_t       y_t


    % state equation
    
    dgp.ABCD.A = zeros(dgp.n_s,dgp.n_s);
    dgp.ABCD.A(1,:) = [dgp.rho, dgp.Sigma * dgp.T^(-dgp.zeta) .* dgp.alpha(2:end)];
    dgp.ABCD.A(3:dgp.n_s,2:dgp.n_s-1) = eye(dgp.n_s-2);
    
    dgp.ABCD.B = [dgp.Sigma * (1 + dgp.T^(-dgp.zeta) * dgp.alpha(1)); 1; zeros(dgp.n_s-2,1)];

    % observation equation
    
    dgp.ABCD.C = dgp.ABCD.A(1,:);
    dgp.ABCD.D = dgp.ABCD.B(1,:);

    % IRFs
    
    dgp.irs_true(i_dgp,:) = compute_irfs(dgp,settings);