% sim_run_setupdgp_arma.m  ARMA: Setup dgp structure for sim_run.m.

indx_rho   = find(dgp_settings.rhos == dgps(1,i_dgp));
indx_theta = find(dgp_settings.thetas == dgps(2,i_dgp));

dgp.n_yp        = dgp.n_y * dgp.p;
dgp.A           = dgp_inputs{indx_rho,indx_theta}.A;
dgp.A_c         = dgp_inputs{indx_rho,indx_theta}.A_c;
dgp.H           = dgp_inputs{indx_rho,indx_theta}.H;
dgp.H_c         = dgp_inputs{indx_rho,indx_theta}.H_c;
dgp.D           = dgp_inputs{indx_rho,indx_theta}.D;
dgp.Sigma       = dgp_inputs{indx_rho,indx_theta}.Sigma;
dgp.alpha_tilde = dgp_inputs{indx_rho,indx_theta}.alpha_tilde;
dgp.n_s         = dgp.n_y * dgp.p + (dgp.n_y * size(dgp.alpha_tilde,1)-dgp.n_eps);