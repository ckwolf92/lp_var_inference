% Shared settings


% Initialize
dgp      = struct;  % Data-generating process (DGP)
settings = struct;  % Estimation settings


% -------------------------------------------------------------------------
% Bootstrap and sample length switches
% -------------------------------------------------------------------------

if boot
    type = ['b' spec];
else
    type = spec;
end

if longT
    type = [type, 'longT'];
    dgp.T = 2000;
else
    dgp.T = 240;  % sample size
end

%--------------------------------------------------------------------------
% dgp and settings setup
%--------------------------------------------------------------------------

% dgp setup
dgp.n_y     = dgp_inputs{1,1}.n_y;
dgp.n_eps   = dgp.n_y;
dgp.T_scale = dgp_settings.T; % artificial sample size to scale mis-specification
dgp.zeta    = dgp_settings.zeta; % mis-specification dampening


% Monte Carlo simulation
settings.sim.numrep      = 5e3; % no. of repetitions
settings.sim.rng_seed    = 202007252; % random number seed
settings.sim.num_workers = 11; % no. of parallel workers (=0: run serial)


% Estimation settings: custom
settings.est.resp_ind  = resp_ind;
settings.est.innov_ind = innov_ind;
settings.est.horzs     = horzs;
settings.est.type      = type; % estimation methods type
settings.est.boot      = boot; % true: bootstrap

% Estimation settings: global
settings.est.no_const  = true; % true: omit intercept
settings.est.se_homosk = true; % true: homoskedastic ses
settings.est.alpha     = 0.1; % significance level
settings.est.boot_num  = 2e3;  % number of bootstrap samples



%--------------------------------------------------------------------------
% specs setup
%--------------------------------------------------------------------------

specs = cell(2,1);           

if settings.est.boot

    specs{1} = {'estimator', 'var',...
            'bias_corr', false,...
            'bootstrap', 'var'};
    specs{2} = {'estimator', 'lp',...
            'bias_corr', false,...
            'bootstrap', 'var'};

else

    specs{1} = {'estimator', 'var',...
            'bias_corr', false};
    specs{2} = {'estimator', 'lp',...
            'bias_corr', false};

end


