function covg = get_asymp_var_covg(dgp, horzs, i, j, alpha)
    
    % get_asymp_var_covg  Get VAR asymptotic coverage
    
    % Inputs:
    % dgp   struct          dgp for simulations
    % horzs [#horizons x 1] Target horizon
    % i     scalar          Response selection index
    % j     scalar          Shock selection index
    
    % Outputs:
    % covg [#horizons x 1]  Coverage probability


avar_VAR  = get_var_avar(dgp, horzs, i, j);
abias_VAR = get_var_abias(dgp, horzs, i, j);
zstar     = norminv(1-alpha/2);
covg      = normcdf(zstar  - abias_VAR./sqrt(avar_VAR)) - ...
    normcdf(-zstar - abias_VAR./sqrt(avar_VAR));


%% Auxiliary functions
    
    function [e_i, e_j] = get_e(dgp, i, j)
        % get_e  Selection vectors for indices i and j

        e_i    = zeros(dgp.n_y*dgp.p, 1);
        e_i(i) = 1;  % Response selection vector
        e_j    = zeros(dgp.n_eps, 1);
        e_j(j) = 1;  % Shock selection vector

    end


    function out = vec(A)
        out = A(:);
    end


    function abias = get_var_abias(dgp, horzs, i, j)
        % get_var_abias  Get VAR asymptotic bias


        % Compute relevant objects
        n_h          = length(horzs);
        [e_i, e_j]   = get_e(dgp, i, j);
        S            = get_S(dgp);
        Psi          = get_psi(dgp, horzs, i, j);
        n_lags_alpha = size(dgp.alpha_tilde, 1)-1;
        abias        = nan(n_h, 1);


        temp1 = zeros(dgp.n_eps, dgp.n_y*dgp.p);  % sum_l alpha_l D H' A'^l-1
        for l=1:n_lags_alpha
            alpha_l = squeeze(dgp.alpha_tilde(l+1,:,:));
            temp1 = temp1 + ...
                alpha_l*dgp.D*dgp.H_c'*(dgp.A_c')^(l-1);
        end

        for jh = 1:n_h  % Loop by horizon

            h     = horzs(jh);
            temp2 = zeros(dgp.n_y*dgp.p, 1);  % sum A^h-l H alpha_l e_j

            if h >=1
                max_l = min(h, n_lags_alpha);
                for l = 1:max_l
                    alpha_l = squeeze(dgp.alpha_tilde(l+1,:,:));
                    temp2   = temp2+dgp.A_c^(h-l)*dgp.H_c*alpha_l*e_j;
                end
            end

            abias(jh) = trace(inv(S)*squeeze(Psi(:,:,jh))* dgp.H_c*temp1) - ...
                e_i'*temp2;
        end

    end


    function avar_VAR = get_var_avar(dgp, horzs, i, j)
        % get_avar_VAR  Get VAR asymptotic variance

        % Preliminaries
        n_h      = length(horzs);
        e_i      = get_e(dgp, i, j);
        Psi      = get_psi(dgp, horzs, i, j);
        S        = get_S(dgp);
        Hbar     = dgp.H_c(:, j+1:end);
        Dbar     = dgp.D(j+1:end, j+1:end);
        avar_VAR = nan(n_h, 1);

        % Loop through each horizon
        for jh = 1:n_h
            h = horzs(jh);
            avar_VAR(jh) = e_i' * dgp.A_c^h * Hbar*Dbar*Hbar'*...
                (dgp.A_c')^h*e_i/(dgp.D(j,j)) +...
                trace(squeeze(Psi(:,:,jh)) * dgp.H_c*dgp.D*...
                dgp.H_c'*squeeze(Psi(:,:,jh))' * inv(S));
        end

    end


    function S = get_S(dgp)
        % get_S  S=Var(ytilde_t)
        % Note: See Assumption 3.1(iv) of the main text.
        
        S = inv(eye(dgp.n_y^2*dgp.p^2) - ...
            kron(dgp.A_c,dgp.A_c))*vec(dgp.H_c*dgp.D*dgp.H_c');
        S = reshape(S, dgp.n_y*dgp.p, []);

    end


    function Psi = get_psi(dgp, horzs, i, j)
        % get_psi  Get Psi_h
        % Note: See Proposition 3.2 of the main text.
    
        n_h = length(horzs);    % Number of horizons
        Psi = zeros(dgp.n_y*dgp.p, dgp.n_y*dgp.p, n_h);
        e_i = get_e(dgp,i,j);  % Response vector

        for jh=1:n_h  % Loop through horizons

            h = horzs(jh);
            if h >= 1
                for l=1:h
                    Psi(:,:,jh) = Psi(:,:,jh) + dgp.A_c^(h-l)*dgp.H_c(:, j) * e_i' * dgp.A_c^(l-1);
                end
            end
        end

    end


end

