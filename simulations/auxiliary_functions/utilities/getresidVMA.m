function [VMA,VARMA] = getresidVMA(VAR_pop,VAR_p,settings)

    % Get Residual VMA Process
    
    % Inputs:
    % VAR_pop   struct          true population VAR(p*)
    % VAR_p     struct          mis-specified VAR(p)
    % settings  struct          estimation settings
    
    % Outputs:
    % VMA       struct          residual VMA process
    % VARMA     struct          implied VARMA(p,\infty)


    %----------------------------------------------------------------
    % Preparations
    %----------------------------------------------------------------
    
    n_y           = size(VAR_pop.Sigma_u,1);
    VMA_hor       = settings.VMA_hor;
    VAR_laglength = settings.VAR_estimlaglength;
    
    %----------------------------------------------------------------
    % VMA Wold Representation
    %----------------------------------------------------------------
    
    IRF_Wold = zeros(n_y,n_y,VMA_hor);
    
    for l = 1:VMA_hor
        for q = 0:min(l-1,VAR_laglength)
            IRF_Wold(:,:,l) = IRF_Wold(:,:,l) + VAR_p.A(:,:,q+1)' * squeeze(VAR_pop.IRF_Wold(l-q,:,:));
        end
    end
    
    VMA.IRF_Wold = permute(IRF_Wold,[3 1 2]);
    
    %----------------------------------------------------------------
    % VARMA Wold Representation
    %----------------------------------------------------------------
    
    % adjust VAR(p) Wold representation
    
    VAR_p.IRF_Wold = zeros(n_y,n_y,VMA_hor);
    VAR_p.IRF_Wold(:,:,1) = eye(n_y);
    
    for l = 1:VMA_hor
    
        if l<VMA_hor
            for j=1:min(l,VAR_laglength)
                VAR_p.IRF_Wold(:,:,l+1) = VAR_p.IRF_Wold(:,:,l+1) + VAR_p.VAR_coeff(1+(j-1)*n_y:j*n_y,:)'*VAR_p.IRF_Wold(:,:,l-j+1);
            end
        end
    
    end
    
    VAR_p.IRF_Wold = permute(VAR_p.IRF_Wold,[3 1 2]);
    
    % get VARMA(p,\infty) Wold representation
    
    VARMA.IRF_Wold = NaN(VMA_hor,n_y,n_y);
    for i_hor = 1:settings.VMA_hor
        VARMA.IRF_Wold(i_hor,:,:) = zeros(1,n_y,n_y);
        for j_hor = 1:i_hor
            VARMA.IRF_Wold(i_hor,:,:) = squeeze(VARMA.IRF_Wold(i_hor,:,:)) + squeeze(VAR_p.IRF_Wold(j_hor,:,:)) * squeeze(VMA.IRF_Wold(i_hor-j_hor+1,:,:));
        end
    end

end