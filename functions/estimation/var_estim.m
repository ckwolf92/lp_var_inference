function [betahat, betahat_varcov, Sigmahat, Sigmahat_varcov, res] ...
            = var_estim(Y, num_lags, se_setting, no_const, df_override)
arguments
Y
num_lags
se_setting
no_const
df_override = []
end
    % Reduced-form Vector Autoregression
    
    % Inputs:
    % Y           T x n      data matrix
    % num_lags    1 x 1      number of lags in VAR
    % se_setting            EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                       OR function handle: function that returns HAC/HAR sandwich matrix
    % no_const    bool      true: omit intercept
    % df_override 1 x 1     Manually override the VAR degree of freedom adjustments (=[] if no override).
    
    % Outputs:
    % betahat           n x (n*num_lags+n+~no_const)      full vector of estimated regression coefficients
    % betahat_varcov    (n x (n*num_lags+n+~no_const))x   var-cov of vec(betahat)
    %                   (n x (n*num_lags+n+~no_const))
    % Sigmahat          n x n                             residual variance-covariance matrix
    % Sigmahat_varcov   (n*(n+1)/2) x (n*(n+1)/2)         var-cov of vech(Sigmahat)
    % res               (T-p) x m                         residuals
    
    % System size
    T = size(Y,1);
    n = size(Y,2);


    % Covariate matrix
    Y_lag = lagmatrix(Y,0:num_lags-1);
    X = Y_lag(num_lags:end-1,:);
    k = size(X,2);
    
    % VAR
    [betahat, betahat_varcov, res, ~] = linreg(Y(num_lags+1:end,:),X,se_setting,no_const);

    % ---------------------------------------------------------------------
    % Degrees of freedom for Sigmahat and Sigmahat_varcov
    % ---------------------------------------------------------------------
    if isempty(df_override)
        df_Sigmahat        = (size(res,1)-n*num_lags-1+no_const);
        df_Sigmahat_varcov = T-(k+(1-no_const));
    else
        df_Sigmahat        = df_override;
        df_Sigmahat_varcov = df_override;
    end
        
    Sigmahat        = (res'*res)/df_Sigmahat;
    D_pl            = D_pl_fn(n);
    Sigmahat_varcov = 2 * D_pl * kron(Sigmahat,Sigmahat) * D_pl'/df_Sigmahat_varcov;

end