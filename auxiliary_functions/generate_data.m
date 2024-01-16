function data_y = generate_data(dgp)

% preparations

T      = dgp.T;

A = dgp.ABCD.A;
B = dgp.ABCD.B;
C = dgp.ABCD.C;
D = dgp.ABCD.D;

n_s   = dgp.n_s;
n_eps = dgp.n_eps;
n_y   = dgp.n_y;

% draw shocks

data_eps = randn(T,n_eps);

% simulate states & measurement error

s = zeros(n_s,1);
data_s = NaN(T,n_s);
for t = 1:T
    s = A * s + B * data_eps(t,:)';
    data_s(t,:) = s';
end

% simulate observables

% data_y = NaN(T, n_y);
% for t = 1:T
%     data_y(t,:) = C * data_s(T_burn + (t-1),:)' + D * data_eps(T_burn + t,:)';
% end

data_y = [0;data_s(:,1)];

end