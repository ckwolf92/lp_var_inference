switch scheme
    case 'mpshock'
        resp_ind     = 2;    % index of response variable of interest
        innov_ind    = 1;    % index of innovation variable
        horzs        = 0:40; % horizons of interest
        SW_model_obs = [10 4]; % (m_shock,y)

    case 'lshock'
        resp_ind     = 2;    
        innov_ind    = 1;    
        horzs        = 0:40; 
        SW_model_obs = [1 19 20 21]; % (shock,pi,w,l)

    case 'mprecursive'
        resp_ind     = 1;    
        innov_ind    = 3;    
        horzs        = 1:40; 
        SW_model_obs = [4 19 5]; % (y,pi,r)

end