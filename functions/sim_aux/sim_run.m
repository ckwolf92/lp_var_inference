% Run simulations

%% PRELIMINARIES

rng(settings.sim.rng_seed, 'twister'); % set RNG seed


numdgp  = size(dgps,2);                 % no. of DGPs
numhorz = length(settings.est.horzs);   % no. of horizons
numspec = length(specs);                % no. of regression specifications
numrep  = settings.sim.numrep;          % no. of repetitions

spec_shared = {'resp_ind',  settings.est.resp_ind, ...
               'innov_ind', settings.est.innov_ind, ...
               'alpha',     settings.est.alpha,...
               'no_const',  settings.est.no_const,...
               'se_homosk', settings.est.se_homosk,...
               'boot_num',  settings.est.boot_num};

dgp.irs_true       = nan(numdgp, numhorz);
dgp.var_asymp_covg = nan(numdgp, numhorz);

%% RUN SIMULATIONS

%----------------------------------------------------------------
% Placeholder Matrices
%----------------------------------------------------------------

estims = zeros(numdgp, numspec, numhorz, numrep);
ses    = estims;

cis_lower = zeros(numdgp, numspec, numhorz, 4, numrep);
cis_upper = cis_lower;

ps = zeros(numdgp, numrep);

%----------------------------------------------------------------
% Start Parallel Workers
%----------------------------------------------------------------

if settings.sim.num_workers > 0
    poolobj = parpool(settings.sim.num_workers);
end

%----------------------------------------------------------------
% Simulations
%----------------------------------------------------------------

for i_dgp = 1:numdgp

    disp(['I am currently at DGP ' num2str(i_dgp)])

    switch exercise
        case 'varma'
            sim_run_setupdgp_sw

        case 'arma'
            sim_run_setupdgp_arma
    end

    dgp.alpha_tilde(2:end,:,:) = dgp.T_scale^(dgp.zeta) * dgp.alpha_tilde(2:end,:,:);
    set_up_varma

    % seed
    i_rand_seeds = randi(2^32-1,1,numrep);

    % Monte Carlo iterations
    parfor (i_rep = 1:numrep, settings.sim.num_workers)

        % seed
        rng(i_rand_seeds(i_rep), 'twister');

        % simulate data

        data_y = generate_data(dgp);

        if settings.est.p_select
            p = ic_var(data_y,settings.p_max,settings.est.p_select);
        else
            p = dgp.p;
        end


        % placeholders
        i_rep_estims    = zeros(numspec, numhorz);
        i_rep_ses       = i_rep_estims;
        i_rep_cis_lower = nan(numspec, numhorz, 4);
        i_rep_cis_upper = i_rep_cis_lower;

        % run specifications

        for i_spec = 1:numspec

            [i_rep_estims(i_spec,:),i_rep_ses(i_spec,:),i_cis_dm,i_cis_boot] ...
                = ir_estim(data_y, p, settings.est.horzs, spec_shared{:}, specs{i_spec}{:});

            i_rep_cis_lower(i_spec,:,1)     = i_cis_dm(1,:);
            i_rep_cis_upper(i_spec,:,1)     = i_cis_dm(2,:);
            i_rep_cis_lower(i_spec,:,2:end) = i_cis_boot(1,:,:);
            i_rep_cis_upper(i_spec,:,2:end) = i_cis_boot(2,:,:);

        end

        % store results

        estims(i_dgp,:,:,i_rep)      = i_rep_estims;
        ses(i_dgp,:,:,i_rep)         = i_rep_ses;
        cis_lower(i_dgp,:,:,:,i_rep) = i_rep_cis_lower;
        cis_upper(i_dgp,:,:,:,i_rep) = i_rep_cis_upper;
        ps(i_dgp,i_rep)              = p;

        if mod(i_rep, ceil(numrep/10)) == 0
            fprintf('%6d%s\n', round(i_rep/numrep*100), '%');
        end

    end

end

clear data_y i_cis_dm i_dgp i_rand_seeds i_rep i_rep_cis_lower i_rep_cis_upper i_rep_estims i_rep_ses ...
    i_spec num2dgp numhorz numrep numspec numdgp spec_shared

if settings.sim.num_workers > 0
    delete(poolobj); % Stop parallel workers
end

%----------------------------------------------------------------
% Store Results
%----------------------------------------------------------------

dgp.dgps = dgps;
clear dgps

results = struct;

results.estims    = estims;
results.ses       = ses;
results.cis_lower = cis_lower;
results.cis_upper = cis_upper;
results.ps        = ps;

clear estims ses cis_lower cis_upper

%----------------------------------------------------------------
% Post-Computations
%----------------------------------------------------------------

% coverage

irs_true_reshape = permute(dgp.irs_true, [1 3 2]);

results.cover_inds = (irs_true_reshape >= results.cis_lower ...
    & irs_true_reshape <= results.cis_upper) + 0; % coverage indicator

results.coverage_prob = mean(results.cover_inds,5); % coverage probability

clear irs_true_reshape

% median length

results.lengths       = results.cis_upper-results.cis_lower;
results.median_length = median(results.lengths,5);

%% SAVE RESULTS

results_filename = ['sim_', exercise, '_', settings.est.type];
save(strcat('results/',results_filename, '.mat'), 'dgp', 'specs', 'settings', 'results');