%% LP/VAR SEs IN RAMEY (2016)
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 07/18/2024

%% HOUSEKEEPING

clear
clc
close all

addpath('auxiliary_functions/')
addpath(genpath('../functions'))

%% ESTIMATION SETTINGS

% shared estimation settings

est.innov_ind = 1;      % shock is ordered first
est.alpha     = 0.1;
est.no_const  = false;  % false: Include intercept
est.se_homosk = true;   % true: homoskedastic ses
est.boot_num  = 2e3;

% application-specific estimation settings

est.resp_ind = [];
est.horzs    = [];
est.p        = [];

% estimator settings

specs{1} = {'estimator', 'var',...
    'bias_corr', false,...
    'bootstrap', 'var'};
specs{2} = {'estimator', 'lp',...
    'bias_corr', false,...
    'bootstrap', 'var'};
numspec = length(specs);

%% SETUP APPLICATIONS

% Gertler-Karadi
appl(1) = setup_appl('Gertler-Karadi', ...
                  readtable('data/Monetarydat.xlsx', 'Sheet', 'Monthly'), ...
                  est, 0:48, 2, @clean_monetary);

% Romer-Romer
appl(2) = setup_appl('Romer-Romer tax', ...
                readtable('data/homtaxdat.xlsx', 'Sheet', 'homtaxdat'), ...
                est, 0:20, 4, @clean_tax);

% Ramey military news
appl(3) = setup_appl('Ramey military news', ...
                readtable('data/homgovdat.xlsx', 'Sheet', 'govdat'), ...
                est, 0:20, 2, @clean_gov);

% FORD TFP
appl(4) = setup_appl('FORD TFP', ...
                readtable('data/Technology_data.xlsx', 'Sheet', 'techdat'),...
                est, 0:20, 2, @clean_tech);

n_appl = length(appl);

%% ESTIMATE IRFs

for i_appl = 1:n_appl  % loop through applications
    
    disp(['Application ' num2str(i_appl) '...'])

    % setup
    appl_i  = appl(i_appl);
    est     = appl_i.est;
    data_y  = [appl_i.data.shock, appl_i.data.y];
    numhorz = length(est.horzs);
    numresp = length(est.resp_ind);
    
    % store results
    estims    = nan(numspec, numresp, numhorz);
    ses       = nan(numspec, numresp, numhorz);
    ses_boot  = nan(numspec, numresp, numhorz);
    cis_lower = nan(numspec, numresp, numhorz, 4);
    cis_upper = nan(numspec, numresp, numhorz, 4);

    tic
    for i_spec = 1:numspec  % LP and VAR loop
        for i_resp = 1:numresp  % outcome variable loop
            
            spec_shared = {'resp_ind', est.resp_ind(i_resp), ...
                'innov_ind', est.innov_ind, ...
                'alpha', est.alpha, 'no_const', est.no_const, ...
                'se_homosk', est.se_homosk, 'boot_num', est.boot_num};

            [irs_i, ses_i, cis_dm_i, cis_boot_i, ses_boot_i] ...
                    = ir_estim(data_y, est.p, est.horzs, spec_shared{:}, specs{i_spec}{:});
            
            % Save
            estims(i_spec, i_resp, :)           = irs_i;
            ses(i_spec, i_resp, :)              = ses_i;
            ses_boot(i_spec, i_resp, :)         = ses_boot_i;
            cis_lower(i_spec, i_resp, :, 1)     = cis_dm_i(1,:);
            cis_lower(i_spec, i_resp, :, 2:end) = cis_boot_i(1,:,:);
            cis_upper(i_spec, i_resp, :, 1)     = cis_dm_i(2,:);
            cis_upper(i_spec, i_resp, :, 2:end) = cis_boot_i(2,:,:);

        end

    end

    runtime = toc;
        
    % save results to appl
    appl(i_appl).results.estims    = estims;
    appl(i_appl).results.ses       = ses;
    appl(i_appl).results.ses_boot  = ses_boot;
    appl(i_appl).results.cis_lower = cis_lower;
    appl(i_appl).results.cis_upper = cis_upper;
    appl(i_appl).results.runtime   = runtime;

end

%% SAVE RESULTS

save('results/res_application.mat', 'appl', '-v7.3')

%% FIGURES

%----------------------------------------------------------------
% IRFs
%----------------------------------------------------------------

% normalizations from Ramey (2016)
scale = {[100, 100, 1, 1]* 0.2567/4.255, ...  
             [1, 1, 1], ...
             [1,1,1]/max(appl(3).results.estims(2, strcmp(appl(3).data.yname, 'g'), :)), ...
             ones(1,3)*100};

for i_appl=1:n_appl
    plot_irfs(appl(i_appl), scale{i_appl})
    exportgraphics(gcf, ['figures/irfs_appl=' num2str(i_appl) '.eps'])
    exportgraphics(gcf, 'figures/irfs_all.pdf', 'Append', i_appl > 1)  % single PDF document
end

%----------------------------------------------------------------
% SEs
%----------------------------------------------------------------

plot_seratio(appl)
exportgraphics(gcf, 'figures/seratio.eps')
exportgraphics(gcf, 'figures/seratio.pdf')  