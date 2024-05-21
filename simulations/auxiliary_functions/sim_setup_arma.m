% ARMA dgp-specific settings

dgp.rhos               = dgp_settings.rhos;
dgp.p                  = dgp_settings.VAR_estimlaglength; % lags
dgp.thetas             = dgp_settings.thetas(1:2); % MA term
settings.est.p_select  = [];

% Set up the DGPs
aux1  = repmat(dgp.rhos,size(dgp.thetas));
aux2  = repmat(dgp.thetas',size(dgp.rhos))';
dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
clear aux1 aux2

