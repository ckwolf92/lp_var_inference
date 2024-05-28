% Smets-Wouters: Setup dgp structure for sim_run.m.

indx_p   = find(dgp_settings.ps == dgps(1,i_dgp));

%----------------------------------------------------------------
% Specification type-specific settings
%----------------------------------------------------------------
switch spec
    case 'fixp'
        settings.est.p  = dgp.ps(i_dgp);
        dgp.alpha_tilde = dgp_inputs{indx_p}.alpha_tilde;

    case 'estp'
        settings.est.p = [];
        dgp.alpha_tilde = dgp_inputs{indx_p}.alpha_tilde;

    case 'worst'
        settings.est.p  = [];
        dgp.alpha_tilde = dgp_inputs{indx_p}.alpha_worst(:,:,:,dgp.worst_hor);

end

%---------------------------------------------------------------------
% Set up DGP
%---------------------------------------------------------------------

dgp.p           = dgp.ps(i_dgp);
dgp.n_yp        = dgp.n_y * dgp.p;
dgp.A           = dgp_inputs{indx_p}.A;
dgp.A_c         = dgp_inputs{indx_p}.A_c;
dgp.H           = dgp_inputs{indx_p}.H;
dgp.H_c         = dgp_inputs{indx_p}.H_c;
dgp.D           = dgp_inputs{indx_p}.D;
dgp.Sigma       = dgp_inputs{indx_p}.Sigma;

dgp.n_s         = dgp.n_y * dgp.p + (dgp.n_y * size(dgp.alpha_tilde,1)-dgp.n_eps);