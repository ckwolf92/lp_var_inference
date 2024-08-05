function VAR_p = popVARp(model,settings,y_aux)

    % Compute VAR(p) in ABCD system
    
    % Inputs:
    % model     struct          ABCD representation of DGP
    % settings  struct          estimation settings
    % y_aux     struct          auxiliary properties of observables y
    
    % Output:
    % VAR_p     struct          implied VAR(p) model

    
    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    n_y = model.n_y;
    
    Sigma_y_big = y_aux.Sigma_y_big;
    
    VAR_laglength   = settings.VAR_estimlaglength;
    VMA_hor         = settings.VMA_hor;
    
    VAR_p.laglength = VAR_laglength;
    
    %----------------------------------------------------------------
    % VAR(p)
    %----------------------------------------------------------------
    
    VAR_p.VAR_coeff = Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,n_y+1:n_y+VAR_laglength*n_y)^(-1) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,1:n_y);
    VAR_p.Sigma_u   = Sigma_y_big(1:n_y,1:n_y) ...
                        - Sigma_y_big(1:n_y,n_y+1:n_y+VAR_laglength*n_y) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,n_y+1:n_y+VAR_laglength*n_y)^(-1) ...
                        * Sigma_y_big(n_y+1:n_y+VAR_laglength*n_y,1:n_y);
    
    VAR_p.A = NaN(n_y,n_y,VAR_laglength+1);
    VAR_p.A(:,:,1) = eye(n_y);
    for l = 2:VAR_laglength+1
        VAR_p.A(:,:,l) = -VAR_p.VAR_coeff(1+(l-2)*n_y:(l-1)*n_y,:);
    end
    
    %----------------------------------------------------------------
    % Wold IRFs
    %----------------------------------------------------------------
    
    VAR_p.IRF_Wold = zeros(n_y,n_y,VMA_hor);
    % VAR_p.IRF_Wold(:,:,1) = VAR_p.Sigma_u^(0.5);
    VAR_p.IRF_Wold(:,:,1) = chol(VAR_p.Sigma_u,'lower');
    
    for l = 1:VMA_hor
        
        if l<VMA_hor
            for j=1:min(l,VAR_laglength)
                VAR_p.IRF_Wold(:,:,l+1) = VAR_p.IRF_Wold(:,:,l+1) + VAR_p.VAR_coeff(1+(j-1)*n_y:j*n_y,:)'*VAR_p.IRF_Wold(:,:,l-j+1);
            end
        end
        
    end
    
    VAR_p.IRF_Wold = permute(VAR_p.IRF_Wold,[3 1 2]);

end