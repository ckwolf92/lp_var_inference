%% GET POPULATION VARMA(p,infty) AS DGP
% this version: 02/07/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

addpath(genpath('../../../functions'))

%% SETTINGS

%----------------------------------------------------------------
% ARMA(1,1) Models
%----------------------------------------------------------------

settings.rho    = 0.9;
settings.thetas = [1/4 1/2 3/4];

settings.n_theta = length(settings.thetas);

%----------------------------------------------------------------
% Computational Settings
%----------------------------------------------------------------

settings.VMA_hor            = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength   = 250; % lag length in population VAR
settings.alpha_lags         = 100;

%----------------------------------------------------------------
% DGP Settings
%----------------------------------------------------------------

settings.ps                 = [1 2 3 4]; % number of lags in mis-specified VAR
settings.max_hor_h          = 10; % maximal IRF horizon of interest
settings.max_hor_alpha_l    = settings.VMA_hor; % maximal lag length for worst-case \alpha(L)
settings.resp_ind           = 1; % response variable of interest
settings.innov_ind          = 1; % innovation variable of interest
settings.T                  = 240; % sample size for DGP
settings.zeta               = 1/2; % mis-specification scaling

settings.n_p                = length(settings.ps);

dgps = cell(settings.n_p,settings.n_theta);

%% GET DGPs

for i_p = 1:settings.n_p

    for i_theta = 1:settings.n_theta

        settings.VAR_estimlaglength = settings.ps(i_p);

        %----------------------------------------------------------------
        % ABCD Representation
        %----------------------------------------------------------------
        
        % parameters
        
        model.param.rho   = settings.rho;
        model.param.theta = settings.thetas(i_theta);
        model.param.sigma = 1;
        
        % ABCD
        
        model.ABCD.A = [model.param.rho, model.param.sigma * model.param.theta; 0, 0];
        model.ABCD.B = [model.param.sigma; 1];
        model.ABCD.C = [model.param.rho, model.param.sigma * model.param.theta];
        model.ABCD.D = model.param.sigma;
        
        % system size
        
        model.n_s   = size(model.ABCD.A,1);
        model.n_eps = size(model.ABCD.B,2);
        model.n_y   = size(model.ABCD.C,1);
        
        %----------------------------------------------------------------
        % VARMA(p,\infty)
        %----------------------------------------------------------------
        
        % VAR(infty)
        
        VAR_infty = popVAR(model,settings);
        y_aux     = get2ndmoments_VAR(VAR_infty,model,settings);
        
        % VAR(p)
        
        VAR_p  = popVARp(model,settings,y_aux);
        
        % Residual VMA(infty)
        
        [VMA,VARMA] = getresidVMA(VAR_infty,VAR_p,model,settings);
        
        %----------------------------------------------------------------
        % Simulation DGP
        %----------------------------------------------------------------
        
        dgps{i_p,i_theta} = dgp_fn(VAR_p,VMA,model,settings);

        clear model VAR_infty VAR_p VARMA VMA y_aux

    end

end

clear i_p i_theta

%% REPORT RESULTS

%----------------------------------------------------------------
% M Table
%----------------------------------------------------------------

% preparations

save_filename = 'm_arma';
theta_names   = {'$\theta = 0.25$','$\theta = 0.50$','$\theta = 0.75$'};

% print table

status = mkdir('tables');
f = fopen(fullfile('tables', strcat(save_filename, '.tex')), 'w'); % open file for writing

fprintf(f, '%s%s%s%s%s\n', '\begin{tabular}{r|', repmat('c', 1, settings.n_theta), '|', repmat('c', 1, settings.n_theta), '}');
fprintf(f, '%s%d%s%d%s\n', '& \multicolumn{', settings.n_theta, '}{c|}{$M$} & \multicolumn{', settings.n_theta, '}{c}{$\frac{M^2}{1+M^2}$} \\');
fprintf(f, '%s', '$p$');
for i=1:2
    for j=1:settings.n_theta
        fprintf(f, '%s%s', ' & ', theta_names{j});
    end
end
fprintf(f, '%s\n%s\n', ' \\', '\hline');

for i_p = 1:settings.n_p
    p = settings.ps(i_p);
    fprintf(f, '%3d', p);
    for i_theta = 1:settings.n_theta
        m = dgps{i_p,i_theta}.M;
        fprintf(f, '%s%5.3f', ' & ', m);
    end
    for i_theta = 1:settings.n_theta
        m2 = dgps{i_p,i_theta}.M2;
        fprintf(f, '%s%5.3f', ' & ', m2);
    end
    fprintf(f, '%s\n', ' \\');
end

fprintf(f, '%s', '\end{tabular}');

fclose(f);

%----------------------------------------------------------------
% \alpha(L) Figure
%----------------------------------------------------------------

status = mkdir('figures');

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [150/255 150/255 150/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];
settings.colors.hor = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255; ... % blue
                            37/255 152/255 14/255]; % green

% lines

settings.line_specs = {'--', ':', '-.'};

% figure selection

h_plot_indx = [1 5 10];
n_h         = length(h_plot_indx);
maxhor_ma   = 20;

indx_p     = [1 4];
indx_theta = [2 2];

% size

plotwidth = 0.41;
gapsize = 0.075;
gapsize_edges = (1-2*plotwidth-gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% figure

figure(1)

for i_p = 1:2

subplot(1,2,i_p)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(i_p);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(1:1:maxhor_ma,dgps{indx_p(i_p),indx_theta(i_p)}.alpha_tilde(2:maxhor_ma+1),'-','Color',settings.colors.black,'LineWidth',4)
for i_h = 1:n_h
    plot(1:1:maxhor_ma,squeeze(dgps{indx_p(i_p),indx_theta(i_p)}.alpha_worst(2:maxhor_ma+1,1,1,h_plot_indx(i_h))),...
        settings.line_specs{i_h},'Color',settings.colors.hor(i_h,:),'LineWidth',4)
end
hold on
xlim([1,maxhor_ma])
% ylim([-1 0.5])
set(gcf,'color','w')
title([' $p$ = ' num2str(settings.ps(indx_p(i_p)))],'interpreter','latex','fontsize',24)
xlabel('$\ell$','interpreter','latex','FontSize',20)
ylabel('$\alpha_\ell$','interpreter','latex','FontSize',20)
if i_p == 2
    legend({'$\alpha(L)$','$\alpha_1^\dagger(L)$', '$\alpha_5^\dagger(L)$','$\alpha_{10}^\dagger(L)$'},...
        'Location','Southeast','fontsize',18,'interpreter','latex','NumColumns',2)
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figures/alpha_arma','-dpng');