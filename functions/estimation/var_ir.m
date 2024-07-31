function [irs, jacob_a, jacob_s] = var_ir_new(A, Sigma, horzs)

    % SVAR(p) impulse responses and Jacobian wrt. parameters
    
    % Inputs:
    % A         n x np  VAR coefficient matrices (A_1, ..., A_p)
    % Sigma     n x n   VAR residual variance-covariance matrix
    % horzs     1 x H   horizons of interest
    
    % Outputs:
    % irs       n x n x H             structural impulse responses Theta_h at select horizons
    % jacob_a   n^2 x (n^2*p) x H     Jacobian of vec(Theta_h) at select horizons wrt. vec(A)
    % jacob_s   n^2 x (n*(n+1)/2) x H Jacobian of vec(Theta_h) at select horizons wrt. vech(Sigma)
    
    
    % Dimensions
    nh     = length(horzs);
    maxh   = max(horzs);
    [n,np] = size(A);
    p      = np/n;

    % Get companion form matrix
    A_c = A_c_fn(A);

    % Structural shock rotation and related auxiliary matrices
    C       = chol(Sigma,'lower');
    C_tilde = C * diag(1./diag(C));  % _tilde: LDL decomposition

    J = [eye(n),zeros(n,n*(p-1))];
    L = L_fn(n);                                          
    K = K_fn(n,n);                                        
    H = L' / (L * (eye(n^2) + K) *  kron_fast(C, eye(n^2),0) * L');  % See e.g. Lutkepohl (2005) Eq 3.7.8
    tmp = diag(-(diag(C).^(-2)));                                    % Jacobian of diag(1./diag(C))
    tmp = tmp(:);  
    H_tilde_aux = diag(tmp) * H;
    H_tilde = kron_fast(C, eye(n^2), 1) * H_tilde_aux + ...
              kron_fast((diag(1./diag(C)))',eye(n^2), 0) * H;  % Product rule   

    % Placeholders
    irs     = zeros(n,n,nh);
    jacob_a = zeros(n^2,n^2*p,nh);
    jacob_s = zeros(n^2,n*(n+1)/2,nh);

    % Lags to aid computations
    ir_p = [C_tilde; zeros(n*(p-1),n)]; % Will contain last p impulse responses, stacked vertically
    jacob_a_p = zeros(n^2,n^2*p,p); % Will contain last p values of the Jacobian of vec(Theta_h) wrt. vec(A)

    % Impact
    if horzs(1) == 0
        irs(:,:,1) = C_tilde;
        Phi_aux = J * J';
        jacob_s(:,:,1) = kron_fast(Phi_aux, eye(n^2), 1) * H_tilde;        
    end

    A_ch = eye(size(A_c));

    % Loop through horizons
    for h=1:maxh
        
        % Compute impulse responses
        A_ch        = A_c*A_ch;
        the_ir      = A_ch*[irs(:,:,1); zeros(n*(p-1),n)];
        the_ir      = the_ir(1:n, 1:n);
        the_past_ir = ir_p(1:n*min(h,p),:);
        ir_p        = [the_ir; ir_p(1:end-n,:)]; % Shift forward in time
        
        the_ind = find(horzs==h,1);
        if ~isempty(the_ind)
            irs(:,:,the_ind) = the_ir; % Store impulse response at select horizons
        end
        
        if nargout>1
            
            % A Jacobian
            the_jacob_a_p                   = zeros(n^2,n^2*p);
            the_jacob_a_p(:,1:n^2*min(h,p)) = kron_fast(the_past_ir', eye(n*n*min(h,p)), 0) ;

            for l=1:min(h,p)
                the_jacob_a_p(:,1:n^2*min(h,p)) = the_jacob_a_p(:,1:n^2*min(h,p)) + kron_fast(A(:,(l-1)*n+1:l*n),jacob_a_p(:,1:n^2*min(h,p),l),1);
            end
            
            jacob_a_p(:,:,2:end) = jacob_a_p(:,:,1:end-1);
            jacob_a_p(:,:,1)     = the_jacob_a_p;
            
            if ~isempty(the_ind)
                jacob_a(:,:,the_ind) = the_jacob_a_p;
            end

            % Sigma Jacobian
            Phi_aux = J * A_ch * J';
            if ~isempty(the_ind)
                jacob_s(:,:,the_ind) = kron_fast(Phi_aux, eye(n^2), 1) * H_tilde;
            end
            
        end
        
    end

end