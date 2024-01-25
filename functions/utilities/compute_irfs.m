function IRF = compute_irfs(dgp, settings)

    % Compute IRFs in ABCD system
    
    % Inputs:
    % A         n_s x n_s       state evolution
    % B         n_s x n_e       state shock response
    % A         n_y x n_s       observables evolution
    % D         n_y x n_e       observables shock response
    % Sigma     n_e x n_e       var-cov of shocks
    % horzs     n_h x 1         IRF horizons of interest
    % innov_ind 1 x 1           index of innovation variable
    % resp_ind  1 x 1           index of response variable
    
    % Output:
    % IRF       n_h x 1         IRFs of response variable to innovation


    % set up system
    
    irf_hor = max(settings.est.horzs+1);
    n_y     = dgp.n_y;
    
    A = dgp.ABCD.A;
    B = dgp.ABCD.B;
    C = dgp.ABCD.C;
    D = dgp.ABCD.D;

    % define shock
    
    if n_y == 1
        shock_weight = 1/dgp.Sigma;
    else
        shock_weight = zeros(n_y,1);
        shock_weight(settings.est.innov_ind,1) = 1/dgp.Sigma(settings.est.innov_ind,settings.est.innov_ind);
    end
    
    % compute IRFs of all observables
    
    IRF = NaN(irf_hor,n_y);
    
    for i = 1:irf_hor
        if i == 1
            IRF(i,:) = (D * shock_weight)';
        else
            IRF(i,:) = C * A^(i-2) * B * shock_weight;
        end
    end

    % select observable of interest
    
    IRF = IRF(settings.est.horzs+1,settings.est.resp_ind);

end