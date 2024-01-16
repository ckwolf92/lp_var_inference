dgp.ABCD.A = zeros(dgp.n_s,dgp.n_s);
dgp.ABCD.A(1,:) = [dgp.rho, dgp.sigma * dgp.T^(-dgp.zeta) .* dgp.alpha(2:end)];
dgp.ABCD.A(3:dgp.n_s,2:dgp.n_s-1) = eye(dgp.n_s-2);

dgp.ABCD.B = [dgp.sigma * (1 + dgp.T^(-dgp.zeta) * dgp.alpha(1)); 1; zeros(dgp.n_s-2,1)];

dgp.ABCD.C = dgp.ABCD.A(1,:);
dgp.ABCD.D = dgp.ABCD.B(1,:);

dgp.irs_true(i_dgp,:) = compute_irfs(dgp,settings);
dgp.irs_true(i_dgp,:) = dgp.irs_true(i_dgp,:) / dgp.sigma;