% sim_setup_sw.m  Smets Wouters-specific setup

switch spec
    case 'fixp'
        dgp.ps                 = dgp_settings.ps;  % lags
        dgp.worst_hor          = [];               % horizon for maximal distortion \alpha(L)
        settings.est.p_select  = [];               % use information criterion to select lag length? (1 for AIC, 2 for BIC)
        settings.p_max         = [];               % maximal lag length to consider

    case 'estp'
        dgp.ps                 = dgp_settings.ps(1);
        dgp.worst_hor          = [];               
        settings.est.p_select  = 1;
        settings.p_max         = 10;

    case 'worst'
        dgp.ps                 = dgp_settings.ps(1);
        dgp.worst_hor          = 4; 
        settings.est.p_select  = 1;
        settings.p_max         = 10;

end

settings.est.ps   = dgp.ps; % lag length used for estimation