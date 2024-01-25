function [irs, ses] = var_select(irs_all, irs_all_varcov, resp_ind, nu)

    % VAR: Return impulse responses of interest along with s.e.
    
    % Inputs:
    % irs_all        n x n x nh      IRFs of all variables to all shocks at all horizons of interest
    % irs_all_varcov n^2 x n^2 x nh  var-cov of IRFs
    % resp_ind       1 x 1           response variable of interest
    % nu             n x 1           shock vector
    
    % Outputs:
    % irs    1 x nh       IRFs of interest
    % ses    1 x nh       st. dev. of IRFs of interest
    

    irs = nu'*permute(irs_all(resp_ind,:,:), [2 3 1]);
    
    % Standard errors
    if nargout>1
        [n,~,nh] = size(irs_all);
        the_eye = eye(n);
        aux = kron(nu',the_eye(resp_ind,:));
        ses = zeros(1,nh);
        for h=1:nh
            ses(h) = sqrt(aux*irs_all_varcov(:,:,h)*aux');
        end
    end

end