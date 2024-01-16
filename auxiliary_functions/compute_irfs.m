function IRF = compute_irfs(dgp,settings);

% prepare

irf_hor = max(settings.est.horzs+1);
n_y     = dgp.n_y;
n_eps   = dgp.n_eps;

A = dgp.ABCD.A;
B = dgp.ABCD.B;
C = dgp.ABCD.C;
D = dgp.ABCD.D;

% compute IRFs

IRF = NaN(irf_hor,n_y);

for i = 1:irf_hor
    if i == 1
        IRF(i,:) = D';
    else
        IRF(i,:) = C * A^(i-2) * B;
    end
end

IRF = IRF(settings.est.horzs+1);

end