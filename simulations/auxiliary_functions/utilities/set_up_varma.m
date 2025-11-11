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
    for i = 1:dgp_settings.alpha_lags
        dgp.ABCD.A(1:dgp.n_yp,dgp.n_yp+1+(i-1)*dgp.n_y:dgp.n_yp+i*dgp.n_y) = dgp.H_c * dgp.T^(-dgp.zeta) * squeeze(dgp.alpha_tilde(i+1,:,:)) * sqrt(dgp.D);
    end
    dgp.ABCD.A(dgp.n_yp+dgp.n_y+1:dgp.n_s,dgp.n_yp+1:dgp.n_s-dgp.n_y) = eye(dgp.n_s-dgp.n_yp-dgp.n_y);
    
    dgp.ABCD.B = zeros(dgp.n_s,dgp.n_eps);
    dgp.ABCD.B(1:dgp.n_yp,:) = dgp.H_c * sqrt(dgp.D);
    dgp.ABCD.B(dgp.n_yp+1:dgp.n_yp+dgp.n_y,:) = eye(dgp.n_y);

    % observation equation
    
    dgp.ABCD.C = dgp.ABCD.A(1:dgp.n_y,:);
    dgp.ABCD.D = dgp.ABCD.B(1:dgp.n_y,:);

    % IRFs
    
    dgp.irs_true(i_dgp,:)        = compute_irfs(dgp,settings);
    dgp.var_asymp_covg(i_dgp, :) = get_asymp_var_covg(dgp, settings.est.horzs, settings.est.resp_ind, ...
        settings.est.innov_ind, settings.est.alpha);