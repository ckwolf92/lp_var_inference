function dgp = dgp_fn(VAR_p,VMA,settings);

    % Get DGP for Simulations
    
    % Inputs:
    % VAR_p     struct          population VAR(p)
    % VMA       struct          residual VMA process
    % settings  struct          estimation settings
    
    % Output:
    % dgp       struct          DGP for simulations

    %----------------------------------------------------------------
    % Parameterization
    %----------------------------------------------------------------

    % system size

    n_y = size(VAR_p.Sigma_u,1);
    
    % autoregressive
    
    dgp.A  = VAR_p.VAR_coeff(1:end-1,:);
    dgp.HD = chol(squeeze(VMA.IRF_Wold(1,:,:)) * squeeze(VMA.IRF_Wold(1,:,:))','lower');
    dgp.H  = zeros(n_y,n_y);
    dgp.D  = zeros(n_y,n_y);
    for i_y = 1:n_y
        dgp.H(:,i_y)   = dgp.HD(:,i_y) / dgp.HD(i_y,i_y);
        dgp.D(i_y,i_y) = dgp.HD(i_y,i_y)^2;
    end
    
    dgp.alpha_tilde = zeros(settings.VMA_hor,n_y,n_y);
    for i_l = 1:settings.VMA_hor
        dgp.alpha_tilde(i_l,:,:) = dgp.H^(-1) * squeeze(VMA.IRF_Wold(i_l,:,:)) * squeeze(VMA.IRF_Wold(1,:,:))^(-1) * dgp.H;
    end
    
    % companion form
    
    dgp.n_y  = size(dgp.A,2);
    dgp.n_yp = dgp.n_y * VAR_p.laglength;
    
    dgp.A_c = A_c_fn(dgp.A');
    dgp.H_c = H_c_fn(dgp.H,VAR_p.laglength);
    
    dgp.Sigma = dgp.H_c * dgp.D * dgp.H_c';
    dgp.S     = reshape((eye(dgp.n_yp^2) - kron(dgp.A_c,dgp.A_c))^(-1) * dgp.Sigma(:),dgp.n_yp,dgp.n_yp);
    
    %----------------------------------------------------------------
    % Properties of \alpha(L)
    %----------------------------------------------------------------
    
    % magnitude of mis-specification
    
    dgp.M_tilde = 0;
    for i_l = 2:settings.VMA_hor
        dgp.M_tilde = dgp.M_tilde + trace(dgp.D * squeeze(dgp.alpha_tilde(i_l,:,:))' * dgp.D^(-1) * squeeze(dgp.alpha_tilde(i_l,:,:)));
    end
    dgp.M_tilde = sqrt(dgp.M_tilde);
    
    clear i_l
    
    dgp.M = settings.T^settings.zeta * dgp.M_tilde;
    dgp.M2 = dgp.M^2/(1 + dgp.M^2);
    
    % worst-case mis-specification
    
    dgp.alpha_worst = alpha_worst_fn(dgp,settings);
    
    % truncate the \alpha(L)'s
    
    dgp.alpha_tilde = dgp.alpha_tilde(1:settings.alpha_lags+1,:,:);
    dgp.alpha_worst = dgp.alpha_worst(1:settings.alpha_lags+1,:,:,:);

end