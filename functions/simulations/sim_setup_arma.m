% sim_setup_arma.m  ARMA dgp-specific settings

dgp.rhos               = dgp_settings.rhos;
dgp.p                  = dgp_settings.VAR_estimlaglength; % lags
dgp.thetas             = dgp_settings.thetas(1:2); % MA term
settings.est.p_select  = [];

