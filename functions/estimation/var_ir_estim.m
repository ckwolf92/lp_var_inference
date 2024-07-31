function [irs, irs_varcov, Ahat_estim, Sigmahat_estim, res_estim] ...
            = var_ir_estim(Y, p, horzs, bias_corr, homosk, no_const, T_eff)
    arguments
        Y
        p
        horzs
        bias_corr
        homosk
        no_const
        T_eff = []
    end
    % VAR(p) least-squares estimates and delta method s.e.
    
    % Inputs:
    % Y           T x n   data vector
    % p           1 x 1   lag length used for impulse response computations
    % horzs       H x 1   horizons of interest
    % bias_corr   bool    true: apply analytical bias correction (Pope, 1990)
    % homosk      bool    true: homoskedastic s.e., false: EHW s.e.
    % no_const    bool    true: omit intercept
    % T_eff       scalar  Specify effective sample size. (=[] for default).

    % Outputs:
    % irs            n x n x H           estimated impulse responses Theta_h at select horizons
    % irs_varcov     n^2 x n^2 x H       var-cov matrices of vec(Theta_h) at select horizons
    % Ahat_estim     n x np              VAR coefficient estimates [A_1,...,A_p] (possibly bias-corrected, possibly including intercept as last column)
    % Sigmahat_estim n x n               VAR residual variance-covariance matrix
    % res_estim      (T-p) x n           estimation residuals
    
    [T,n] = size(Y);
    
    % One-step forecasting regression of Y_{t+1} on (Y_t, ..., Y_{t-p+1})
    [Ahat_estim, Ahat_estim_varcov, Sigmahat_estim, Sigmahat_estim_varcov, res_estim] ...
        = var_estim(Y, p, homosk, no_const, T_eff);
    
    % Bias correction, if desired
    if bias_corr
        Ahat_estim(:,1:end-1+no_const) = var_biascorr(Ahat_estim(:,1:end-1+no_const), Sigmahat_estim, T);
    end
    
    % Only use first p VAR coefficient matrices to compute impulse responses
    Ahat = Ahat_estim(:,1:n*p);
    Ahat_varcov = Ahat_estim_varcov(1:n^2*p,1:n^2*p);
    Sigmahat = Sigmahat_estim;
    Sigmahat_varcov = Sigmahat_estim_varcov;
    
    if nargout==1
        irs = var_ir(Ahat,Sigmahat,horzs); % Compute impulse responses
    else
        [irs, jacob_a, jacob_s] = var_ir(Ahat,Sigmahat,horzs); % Compute impulse responses and Jacobian
        nh = length(horzs);
        irs_varcov = zeros(n^2,n^2,nh);
        for h=1:nh
            irs_varcov(:,:,h) = jacob_a(:,:,h)*Ahat_varcov*jacob_a(:,:,h)' ...
                + jacob_s(:,:,h)*Sigmahat_varcov*jacob_s(:,:,h)';
        end
    end

end