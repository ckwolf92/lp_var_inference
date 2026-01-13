function VAR_p = var_fn(data,laglength,settings)

    % Compute VAR(p) given dataset
    
    % Inputs:
    % data      T x n      data matrix
    % laglength 1 x 1      lag length
    % settings  struct     settings
    
    % Output:
    % VAR_p     struct          implied VAR(p) model

    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------

    VAR_p.laglength   = laglength;
    VAR_laglength     = VAR_p.laglength;
    VMA_hor           = settings.VMA_hor;
    n_y               = size(data,2);

    %----------------------------------------------------------------
    % Estimate VAR
    %----------------------------------------------------------------

    [VAR_p.VAR_coeff, ~, VAR_p.Sigma_u, ~] ...
            = var_estim(data, VAR_p.laglength, false, false);

    VAR_p.VAR_coeff = VAR_p.VAR_coeff';

    VAR_p.A = NaN(n_y,n_y,VAR_laglength+1);
    VAR_p.A(:,:,1) = eye(n_y);
    for l = 2:VAR_laglength+1
        VAR_p.A(:,:,l) = -VAR_p.VAR_coeff(1+(l-2)*n_y:(l-1)*n_y,:);
    end

    %----------------------------------------------------------------
    % Wold IRFs
    %----------------------------------------------------------------
    
    VAR_p.IRF_Wold = zeros(n_y,n_y,VMA_hor);
    VAR_p.IRF_Wold(:,:,1) = chol(VAR_p.Sigma_u,'lower');
    
    for l = 1:VMA_hor
        
        if l<VMA_hor
            for j=1:min(l,VAR_laglength)
                VAR_p.IRF_Wold(:,:,l+1) = VAR_p.IRF_Wold(:,:,l+1) + VAR_p.VAR_coeff(1+(j-1)*n_y:j*n_y,:)'*VAR_p.IRF_Wold(:,:,l-j+1);
            end
        end
        
    end
    
    VAR_p.IRF_Wold = permute(VAR_p.IRF_Wold,[3 1 2]);