function alpha_worst = alpha_worst_fn(dgp,settings);

    % Get DGP for Simulations
    
    % Inputs:
    % dgp         struct          simulation DGP
    % settings    struct          estimation settings
    
    % Output:
    % alpha_worst struct          worst-case \alpha(L)

    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    n_y  = dgp.n_y;
    n_yp = dgp.n_yp;
    
    max_hor_alpha_l = settings.max_hor_alpha_l;
    max_hor_h       = settings.max_hor_h;
    h_list          = [1:1:max_hor_h];
    
    resp_ind  = settings.resp_ind;
    innov_ind = settings.innov_ind;
    
    e_i   = zeros(n_yp,1); e_i(resp_ind,1) = 1;
    e_j   = zeros(n_y,1);  e_j(innov_ind,1) = 1;
    
    A = dgp.A_c;
    D = dgp.D;
    H = dgp.H_c;
    S = dgp.S;
    
    M_tilde     = dgp.M_tilde;
    alpha_tilde = dgp.alpha_tilde;
    
    %----------------------------------------------------------------
    % Auxiliary Objects
    %----------------------------------------------------------------
    
    Psi = NaN(n_yp,n_yp,max_hor_h);
    for h = 1:max_hor_h
        Psi(:,:,h) = zeros(n_yp,n_yp);
        for l = 1:h
            Psi(:,:,h) = Psi(:,:,h) + A^(h-l) * H(:,innov_ind) * e_i' * A^(l-1);
        end
    end
    
    %----------------------------------------------------------------
    % Xi
    %----------------------------------------------------------------
    
    xi = NaN(n_y,n_y,max_hor_alpha_l,max_hor_h);
    
    for i = 1:max_hor_h
        for j = 1:max_hor_alpha_l
            h = h_list(i);
            l = j;
            if l <= h
                xi(:,:,j,i) = D^(0.5) * H' * (A')^(l-1) * S^(-1) * Psi(:,:,h) * H * D^(0.5) ...
                                - (D^(-0.5) * e_j * e_i' * A^(h-l) * H * D^(0.5));
            else
                xi(:,:,j,i) = D^(0.5) * H' * (A')^(l-1) * S^(-1) * Psi(:,:,h) * H * D^(0.5);
            end
        end
    end
    
    %----------------------------------------------------------------
    % \alpha(L)
    %----------------------------------------------------------------
    
    alpha_worst = zeros(n_y,n_y,max_hor_alpha_l,max_hor_h);
    
    for i = 1:max_hor_h
        for j = 1:max_hor_alpha_l-1
            alpha_worst(:,:,j+1,i) = xi(:,:,j,i)';
        end
    end
    
    for i = 1:max_hor_h
        alpha_worst_norm = 0;
        for i_l = 1:max_hor_alpha_l
            alpha_worst_norm = alpha_worst_norm + trace(D * squeeze(alpha_worst(:,:,i_l,i))' ...
                * D^(-1) * squeeze(alpha_worst(:,:,i_l,i)));
        end
        alpha_worst_norm = sqrt(alpha_worst_norm);
        scale_i = M_tilde/alpha_worst_norm;
        alpha_worst(:,:,:,i) = -alpha_worst(:,:,:,i) * scale_i;
    end
    
    alpha_worst = permute(alpha_worst,[3 1 2 4]);

end