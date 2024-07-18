%% WORST-CASE FIGURES
% Jose L. Montiel Olea, Mikkel Plagborg-Moller, Eric Qian, and Christian Wolf
% this version: 07/18/2024

%% HOUSEKEEPING

clear
clc
close all

addpath(genpath('../functions'))
addpath(genpath('../emp_ses'))

%% SETTINGS

%----------------------------------------------------------------
% Worst-Case Computations
%----------------------------------------------------------------

signif = 0.1;                           % significance level
se_ratio = linspace(0.001,0.999,100);   % values of SE ratio to plot (sqrt{aVar(AR)/aVar(LP)})
M = [0.1 1 1.5 2 3];                    % values of M to plot
linestyles = {'-','--','-.',':','.'};   % line style for each value of M
plot_pos = [0 0 8 4];                   % figure position (in inches)
plot_ylim_length = [0 2];               % vertical axis limits for relative length plot

%----------------------------------------------------------------
% Empirical SEs
%----------------------------------------------------------------

% import results

load res_application

% select variables and horizons

var_select = {[1 2 3 4]; ...
                [1 2 3]; ...
                [1 2 3]; ...
                [1 2 3]; ...
                [1 2 3 4]; ...
                [1 2 3]; ...
                [1 2 3]; ...
                [1 2 3]};

hor_select = {[1:1:49]; ...
                [1:1:21]; ...
                [1:1:21]; ...
                [1:1:21];...
                [13:1:49]; ...
                [5:1:21]; ...
                [5:1:21]; ...
                [5:1:21]};

n_appl = length(var_select);

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.3 * settings.colors.blue + 0.7 * [1 1 1];

%% WORST-CASE COMPUTATIONS

%----------------------------------------------------------------
% Set-Up
%----------------------------------------------------------------

n_M = length(M);
n_r = length(se_ratio);

z = norminv(1-signif/2); % normal critical value
r = @(b,c) normcdf(b-c)+normcdf(-b-c); % rejection probability
cv = @(b) fzero(@(c) r(b,c)-signif,z+b); % bias-aware critical value
prob_fct = @(k,b) log(r(b,z))+log(1-r(b/k,z)); % log joint probability Hausman rejects and VAR CI doesn't cover

%----------------------------------------------------------------
% Get Worst-Case Results
%----------------------------------------------------------------

% coverage

wc_cov = 1-r(M(:).*sqrt(1./(se_ratio(:).^2)'-1),z);

% bias-aware + Hausman test

cv_var_aware = nan(n_M,n_r);
weight_opt_aware = nan(n_M,n_r);
length_opt_aware = nan(n_M,n_r);
wc_prob = nan(1,n_r);
wc_bias = nan(1,n_r);
opts = optimoptions('fmincon','Display','notify');

for i_r = 1:n_r

    the_aux = sqrt(1./se_ratio(i_r)^2-1);

    for i_M=1:n_M

        % critical value for bias-aware VAR CI
        cv_var_aware(i_M,i_r) = cv(M(i_M)*the_aux);

        % optimal bias-aware CI
        the_obj = @(w) cv((1-w)*M(i_M)*the_aux/sqrt(1+w^2*the_aux^2))*sqrt(1+w^2*the_aux^2);
        [the_w,the_val] = fmincon(the_obj, M(i_M)^2/(1+M(i_M)^2), [], [], [], [], 0, 1, [], opts);
        weight_opt_aware(i_M,i_r) = the_w;
        length_opt_aware(i_M,i_r) = the_val*se_ratio(i_r)/z;

    end

    % worst-case asy. probability of the joint event:
    % (a) Hausman test fails to reject correct VAR specification
    % AND (b) conventional VAR t-test rejects true parameter
    [the_b,the_val] = fmincon(@(b) -prob_fct(the_aux,b), the_aux, [], [], [], [], 1e-9, Inf, [], opts);
    wc_prob(i_r) = exp(-the_val);
    wc_bias(i_r) = the_b;

end

length_var_aware = (cv_var_aware/z).*se_ratio; % relative length of bias-aware VAR CI

%% EXTRACT SE RATIOS

%----------------------------------------------------------------
% Ratios
%----------------------------------------------------------------

ses     = cell(n_appl,1);

for i_appl = 1:n_appl

    ses_tmp = NaN(length(hor_select{i_appl}),length(var_select{i_appl}));
    for i_hor = 1:length(hor_select{i_appl})
        for i_var = 1:length(var_select{i_appl})
            if i_appl <= 4
                ses_tmp(i_hor,i_var) = appl(i_appl).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    ./ appl(i_appl).results.ses_boot(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
            else
                ses_tmp(i_hor,i_var) = appl(i_appl-4).results.ses_boot(1,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor)) ...
                    ./ appl(i_appl-4).results.ses_boot(2,var_select{i_appl}(i_var),hor_select{i_appl}(i_hor));
            end
        end
    end
    ses_tmp = ses_tmp(:);
    ses{i_appl} = ses_tmp;
    clear ses_tmp

end

%----------------------------------------------------------------
% Percentiles
%----------------------------------------------------------------

ses_all_1 = [ses{1};ses{2};ses{3};ses{4}]; % all horizons
ses_all_2 = [ses{5};ses{6};ses{7};ses{8}]; % medium/long horizons

ses_lb_1 = prctile(ses_all_1,10);
ses_ub_1 = prctile(ses_all_1,90);
ses_lb_2 = prctile(ses_all_2,10);
ses_ub_2 = prctile(ses_all_2,90);

%% PLOT RESULTS

%----------------------------------------------------------------
% Worst-Case Coverage of VAR CI
%----------------------------------------------------------------

for i_figure = 1:2

    figure('Units','inches','Position',plot_pos);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    for i_M=1:n_M
        plot(se_ratio,wc_cov(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
        hold on
    end
    if i_figure > 1
        jbfill([ses_lb_2 ses_ub_2],[0 0],[1 1],...
            settings.colors.lblue,settings.colors.lblue,0,0.5);
        hold on
    end
    the_xlim = [0 1];
    line(the_xlim,[1 1]*(1-signif),'Color','k','LineStyle','-');
    hold on
    xlim(the_xlim);
    ylim([0 1]);
    yticks([0:0.2:1])
    label_legend(M);
    print(['results/wc_coverage_' num2str(i_figure)],'-dpng');

end

%----------------------------------------------------------------
% Relative Length of Bias-Aware VAR CI
%----------------------------------------------------------------

for i_figure = 1:2

    figure('Units','inches','Position',plot_pos);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    for i_M=1:n_M
        plot(se_ratio,length_var_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
        hold on
    end
    if i_figure > 1
        jbfill([ses_lb_2 ses_ub_2],[0 0],[2 2],...
            settings.colors.lblue,settings.colors.lblue,0,0.5);
        hold on
    end
    the_xlim = [0 1];
    line(the_xlim,[1 1],'Color','k','LineStyle','-');
    hold on
    xlim(the_xlim);
    ylim([0 2]);
    yticks([0:0.5:2])
    label_legend(M);
    print(['results/ba_rellength_' num2str(i_figure)],'-dpng');

end

%----------------------------------------------------------------
% LP Weight in Optimal Bias-Aware CI
%----------------------------------------------------------------

for i_figure = 1:2

    figure('Units','inches','Position',plot_pos);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    for i_M=1:n_M
        plot(se_ratio,weight_opt_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
        hold on
    end
    if i_figure > 1
        jbfill([ses_lb_2 ses_ub_2],[0 0],[1 1],...
            settings.colors.lblue,settings.colors.lblue,0,0.5);
        hold on
    end
    the_xlim = [0 1];
    xlim(the_xlim);
    ylim([0 1]);
    yticks([0:0.2:1])
    label_legend(M);
    print(['results/ba_lpoptweight_' num2str(i_figure)],'-dpng');

end

%----------------------------------------------------------------
% Relative Length of Optimal Bias-Aware CI
%----------------------------------------------------------------

for i_figure = 1:2

    figure('Units','inches','Position',plot_pos);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    for i_M=1:n_M
        plot(se_ratio,length_opt_aware(i_M,:),linestyles{i_M},'Color','k','LineWidth',1);
        hold on
    end
    if i_figure > 1
        jbfill([ses_lb_2 ses_ub_2],[0 0],[1 1],...
            settings.colors.lblue,settings.colors.lblue,0,0.5);
        hold on
    end
    the_xlim = [0 1];
    xlim(the_xlim);
    ylim([0 1]);
    yticks([0:0.2:1])
    label_legend(M);
    print(['results/ba_optrellength_' num2str(i_figure)],'-dpng');

end

%----------------------------------------------------------------
% Worst-Case Joint Probability
%----------------------------------------------------------------

for i_figure = 1:2

    figure('Units','inches','Position',plot_pos);
    set(gca,'TickLabelInterpreter','latex')
    hold on
    if i_figure > 1
        jbfill([ses_lb_2 ses_ub_2],[0 0],[1 1],...
            settings.colors.lblue,settings.colors.lblue,0,0.5);
        hold on
    end
    plot(se_ratio,wc_prob,'Color','k','LineWidth',2);
    the_xlim = [0 1];
    line(the_xlim,[1 1]*signif,'Color','k','LineStyle',':');
    hold on
    xlim(the_xlim);
    ylim([0 1]);
    yticks([0:0.2:1])
    label_legend();
    print(['results/wc_joint_' num2str(i_figure)],'-dpng');

end

%% PLOTTING FUNCTION

function label_legend(varargin)
    xlabel('$\sqrt{\mathrm{aVar}(\hat{\delta}_h)/\mathrm{aVar}(\hat{\beta}_h)}$','Interpreter','Latex');
    set(gca,'FontSize',12);
    if ~isempty(varargin)
        legend(strcat('M=', arrayfun(@num2str, varargin{1}, 'UniformOutput', 0)),...
            'Interpreter','latex','Location','SouthEast','FontSize',12);
    end
end