% Bring VAR model into state-space form & compute IRFs

    % Original system, VAR(1) companion form:
    % Y_t = A_c Y_{t-1} + (H + T^{-\zeta} * \alpha(L)) * \epsilon_t
    % where
    % Y         (y_t, y_{t-1}, ...)
    
    % State-space form
    % s_t = A s_{t-1} + B \epsilon_t
    % x_t = C s_{t-1} + D \epsilon_t
    % where
    % s_t       (Y_t, \epsilon_t, \epsilon_{t-1}, ...)'
    % x_t       y_t


    % state equation
    
    dgp.ABCD.A = zeros(dgp.n_s,dgp.n_s);
    dgp.ABCD.A(1:dgp.n_yp,1:dgp.n_yp) = dgp.A_c;
    dgp.ABCD.A(1:dgp.n_y,dgp.n_yp+1:end) = (dgp.T^(-dgp.zeta) .* dgp.alpha(:,dgp.n_y+1:end)) * blkdiag(dgp.Sigma_aux{:});
    dgp.ABCD.A(dgp.n_yp+dgp.n_y+1:dgp.n_s,dgp.n_yp+1:dgp.n_s-dgp.n_y) = eye(dgp.n_s-dgp.n_yp-dgp.n_y);
    
    dgp.ABCD.B = zeros(dgp.n_s,dgp.n_eps);
    dgp.ABCD.B(1:dgp.n_y,:) = dgp.H * dgp.Sigma;
    dgp.ABCD.B(dgp.n_yp+1:dgp.n_yp+dgp.n_y,:) = eye(dgp.n_y);

    % observation equation
    
    dgp.ABCD.C = dgp.ABCD.A(1:dgp.n_y,:);
    dgp.ABCD.D = dgp.ABCD.B(1:dgp.n_y,:);

    % IRFs
    
    dgp.irs_true(i_dgp,:) = compute_irfs(dgp,settings);