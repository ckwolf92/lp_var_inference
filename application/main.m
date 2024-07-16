%% PRELIMINARIES

clear
clc
close all

addpath('auxiliary_functions/')
addpath(genpath('../functions'))

%% ESTIMATION SETTINGS

% Shared estimation settings
est.innov_ind = 1;      % Shock is ordered first
est.alpha     = 0.1;
est.no_const  = false;  % false: Include intercept
est.se_homosk = true;   % true: homoskedastic ses
est.boot_num  = 1e1;

% Application-specific estimation settings
est.resp_ind = [];
est.horzs    = [];
est.p        = [];

% Estimator settings
specs{1} = {'estimator', 'var',...
    'bias_corr', false,...
    'bootstrap', 'var'};
specs{2} = {'estimator', 'lp',...
    'bias_corr', false,...
    'bootstrap', 'var'};
numspec = length(specs);

%% SETUP APPLICATIONS

% Gertler-Karadi
Appl(1) = setup_Appl('Gertler-Karadi', ...
                  readtable('Ramey_HOM_data/Monetarydat.xlsx', 'Sheet', 'Monthly'), ...
                  est, 0:48, 2, @clean_monetary);

% Romer-Romer
Appl(2) = setup_Appl('Romer-Romer tax', ...
    readtable('Ramey_HOM_data/homtaxdat.xlsx', 'Sheet', 'homtaxdat'), ...
    est, 0:20, 4, @clean_tax);

% Ramey military news
Appl(3) = setup_Appl('Ramey military news', ...
readtable('Ramey_HOM_data/homgovdat.xlsx', 'Sheet', 'govdat'), ...
    est, 0:20, 2, @clean_gov);

% FORD TFP
Appl(4) = setup_Appl('FORD TFP', ...
    readtable('Ramey_HOM_data/Technology_data.xlsx', 'Sheet', 'techdat'),...
    est, 0:20, 2, @clean_tech);


numappl = length(Appl);

%% ESTIMATE IRFs


for i_appl = 1:numappl  % Loop through applications
    
    disp(['Application ' num2str(i_appl) '...'])

    % Setup
    Appl_i  = Appl(i_appl);
    est     = Appl_i.est;
    data_y  = [Appl_i.data.shock, Appl_i.data.y];
    numhorz = length(est.horzs);
    numresp = length(est.resp_ind);
    
    % Store results
    estims    = nan(numspec, numresp, numhorz);
    ses       = nan(numspec, numresp, numhorz);
    cis_lower = nan(numspec, numresp, numhorz, 4);
    cis_upper = nan(numspec, numresp, numhorz, 4);

    tic
    for i_spec = 1:numspec  % LP and VAR loop
        for i_resp = 1:numresp  % Outcome variable loop
            
            spec_shared = {'resp_ind', est.resp_ind(i_resp), ...
                'innov_ind', est.innov_ind, ...
                'alpha', est.alpha, 'no_const', est.no_const, ...
                'se_homosk', est.se_homosk, 'boot_num', est.boot_num};


            [irs_i, ses_i, cis_dm_i, cis_boot_i] ...
                = ir_estim(data_y, est.p, est.horzs, spec_shared{:}, specs{i_spec}{:});
            
            % Save
            estims(i_spec, i_resp, :)           = irs_i;
            ses(i_spec, i_resp, :)              = ses_i;
            cis_lower(i_spec, i_resp, :, 1)     = cis_dm_i(1,:);
            cis_lower(i_spec, i_resp, :, 2:end) = cis_boot_i(1,:,:);
            cis_upper(i_spec, i_resp, :, 1)     = cis_dm_i(2,:);
            cis_upper(i_spec, i_resp, :, 2:end) = cis_boot_i(2,:,:);


        end

    end

    runtime = toc;
        
    % Save results to Appl
    Appl(i_appl).results.estims    = estims   ;
    Appl(i_appl).results.ses       = ses      ;
    Appl(i_appl).results.cis_lower = cis_lower;
    Appl(i_appl).results.cis_upper = cis_upper;
    Appl(i_appl).results.runtime   = runtime;

end

save('results/res_application.mat', 'Appl', '-v7.3')

%% Plot IRFs

% IRF normalizations from Ramey (2016)
scale = {[100, 100, 1, 1]* 0.2567/4.255, ...  
             [1, 1, 1], ...
             [1,1,1]/max(Appl(3).results.estims(2, strcmp(Appl(3).data.yname, 'g'), :)), ...
             ones(1,3)*100};


for i_appl=1:numappl
    plot_irfs(Appl(i_appl), scale{i_appl})
    exportgraphics(gcf, ['figures/irfs_Appl=' num2str(i_appl) '.eps'])
    exportgraphics(gcf, 'figures/irfs_all.pdf', 'Append', i_appl > 1)  % Single PDF document
end

%% Plot SEs

plot_seratio(Appl)
exportgraphics(gcf, 'figures/seratio.eps')
exportgraphics(gcf, 'figures/seratio.pdf')  

function plot_seratio(Appl)
%%
    numappl = length(Appl);
    f = figure('Visible', 'off', 'Units','inches');
    f.Position(3:4) = [6,4];
    ax              = gca;
    ax.FontSize     = 12;

    Style = {'-', '--', ':', '-.'};

    for i_appl = 1:numappl
        se_ratio = squeeze(Appl(i_appl).results.ses(1,:,:) ./ ...
                           Appl(i_appl).results.ses(2,:,:))';
        horzs = Appl(i_appl).est.horzs;
        subplot(2,2,i_appl)
        plot(horzs([1,end]), [1,1], 'Color', 'k', 'LineWidth', 1, ...
            'LineStyle', '-')
        hold on
        
        p = [];


        for i_plot = 1:size(se_ratio, 2)
            p_i = plot(horzs, se_ratio(:, i_plot), 'LineWidth', 1, ...
                'Color', 'b','LineStyle', Style{i_plot});
            p = [p; p_i];
            ylim([0, 1.1*max(se_ratio(:))])
        end
        

        % Visuals
        legend(p, Appl(i_appl).data.yname, 'Location', 'southoutside', ...
            'Interpreter','latex', 'Orientation','horizontal', ...
            'FontSize', 9, 'NumColumns',2)
        title(Appl(i_appl).name, 'Interpreter','latex')
        box on; grid on;
        xlim([horzs(1) horzs(end)])

        % Set ticks
        ticksize = 12/(month(Appl(i_appl).data.date(2)) - ...
                       month(Appl(i_appl).data.date(1)));
        xticks(horzs(1):ticksize:horzs(end))
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        sgtitle('SE ratio: VAR to LP', 'Interpreter', 'latex')
    end


end