clear all;

% Worst-case coverage of conventional VAR CI/ellipsoid,
% comparison of bias-aware VAR CI to LP CI,
% and optimal bias-aware CI


%% Settings

signif = 0.1;                           % Significance level
se_ratio = sort([linspace(0.001,0.999,100) 0.25 0.5]); % Values of SE ratio to plot (sqrt{aVar(AR)/aVar(LP)})
k = [1 2 5 10 20];                      % Values of k to plot (for joint Wald confidence ellipsoid)
M = [0.1 1 1.5 2 3];                    % Values of M to plot
linestyles = {'-k','--k','-.k',':k','.k'}; % Line style for each value of M
plot_pos = [0 0 8 4];                   % Figure position (in inches)
plot_ylim_length = [0 2];               % Vertical axis limits for relative length plot
plot_ellipse_k = [2 3 4];               % Which indices of the k vector to use as subplots when plotting worst-case ellipsoid coverage as a function of M
plot_ellipse_M = [2 3 4];               % Which indices of the M vector to use as subplots when plotting worst-case ellipsoid coverage as a function of k


%% Compute

n_M = length(M);
n_r = length(se_ratio);

z = norminv(1-signif/2); % Normal critical value
r = @(b,c) normcdf(b-c)+normcdf(-b-c); % Rejection probability
cv = @(b) fzero(@(c) r(b,c)-signif,z+b); % Bias-aware critical value (note: this is faster than ncx2inv command)
prob_fct = @(k,b) log(r(b,z))+log(1-r(b/k,z)); % Log joint probability Hausman rejects and VAR CI doesn't cover

noncentr = M(:).^2.*(se_ratio(:).^(-2)'-1); % Worst-case noncentrality parameter
wc_cov = 1-r(sqrt(noncentr),z); % Worst-case coverage of conventional CI as function of relative length

cv_var_aware = nan(n_M,n_r);
weight_opt_aware = nan(n_M,n_r);
length_opt_aware = nan(n_M,n_r);
wc_prob = nan(1,n_r);
wc_bias = nan(1,n_r);
opts = optimoptions('fmincon','Display','notify');

for i_r=1:n_r

    the_aux = sqrt(1./se_ratio(i_r)^2-1);

    for i_M=1:n_M

        % Critical value for bias-aware VAR CI
        cv_var_aware(i_M,i_r) = cv(sqrt(noncentr(i_M,i_r)));

        % Optimal bias-aware CI
        the_obj = @(w) cv((1-w)*M(i_M)*the_aux/sqrt(1+w^2*the_aux^2))*sqrt(1+w^2*the_aux^2);
        [the_w,the_val] = fmincon(the_obj, M(i_M)^2/(1+M(i_M)^2), [], [], [], [], 0, 1, [], opts);
        weight_opt_aware(i_M,i_r) = the_w;
        length_opt_aware(i_M,i_r) = the_val*se_ratio(i_r)/z;

    end

    % Worst-case asy. probability of the joint event:
    % (a) Hausman test fails to reject correct VAR specification
    % AND (b) conventional VAR t-test rejects true parameter
    [the_b,the_val] = fmincon(@(b) -prob_fct(the_aux,b), the_aux, [], [], [], [], 1e-9, Inf, [], opts);
    wc_prob(i_r) = exp(-the_val);
    wc_bias(i_r) = the_b;

end

length_var_aware = (cv_var_aware/z).*se_ratio; % Relative length of bias-aware VAR CI

% Worst-case coverage of joint Wald confidence ellipsoid
n_k = length(k);
wc_prob_wald = nan(n_M,n_r,n_k);
for i_k=1:n_k
    wc_prob_wald(:,:,i_k) = ncx2cdf(chi2inv(1-signif,k(i_k)),k(i_k),noncentr);
end


%% Plot

% Worst-case coverage of VAR CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,wc_cov(i_M,:),linestyles{i_M},'LineWidth',1);
end
hold off;
the_xlim = xlim;
line(the_xlim,[1 1]*(1-signif),'Color','k','LineStyle','-'); % Mark nominal coverage level on vertical axis
xlim(the_xlim);
ylim([0 1]);
label_axis(true);
label_legend(M,'M');

% Relative length of bias-aware VAR CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,length_var_aware(i_M,:),linestyles{i_M},'LineWidth',1);
end
hold off;
the_xlim = xlim;
line(the_xlim,[1 1],'Color','k','LineStyle','-'); % Mark 1 on vertical axis
xlim(the_xlim);
ylim(plot_ylim_length);
label_axis(true);
label_legend(M,'M');

% Weight on LP in optimal bias-aware CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,weight_opt_aware(i_M,:),linestyles{i_M},'LineWidth',1);
end
hold off;
ylim([0 1]);
label_axis(true);
label_legend(M,'M');

% Relative length of optimal bias-aware CI
figure('Units','inches','Position',plot_pos);
hold on;
for i_M=1:n_M
    plot(se_ratio,length_opt_aware(i_M,:),linestyles{i_M},'LineWidth',1);
end
hold off;
ylim([0 1]);
label_axis(true);
label_legend(M,'M');

% Worst-case joint probability
figure('Units','inches','Position',plot_pos);
plot(se_ratio,wc_prob,'Color','k','LineWidth',2);
the_xlim = xlim;
line(the_xlim,[1 1]*signif,'Color','k','LineStyle',':'); % Mark nominal significance level on vertical axis
xlim(the_xlim);
ylim([0 1]);
label_axis(true);

% Worst-case coverage of Wald confidence ellipsoid
n_plot_k = length(plot_ellipse_k);
figure('Units','inches','Position',plot_pos.*[1 1 1.5 1]);
for i_k=1:n_plot_k
    subplot(1,n_plot_k,i_k);
    hold on;
    for i_M=1:n_M
        plot(se_ratio,wc_prob_wald(i_M,:,plot_ellipse_k(i_k)),linestyles{i_M},'LineWidth',1);
    end
    hold off;
    the_xlim = xlim;
    line(the_xlim,[1 1]*(1-signif),'Color','k','LineStyle','-'); % Mark nominal coverage level on vertical axis
    xlim(the_xlim);
    ylim([0 1]);
    label_axis(false);
    title(sprintf('%s%d%s','$k=',k(plot_ellipse_k(i_k)),'$'),'Interpreter','Latex','FontSize',14);
    if i_k==n_plot_k
        label_legend(M,'M');
    end
end

figure('Units','inches','Position',plot_pos.*[1 1 1.5 1]);
n_plot_M = length(plot_ellipse_M);
for i_M=1:n_plot_M
    subplot(1,n_plot_M,i_M);
    hold on;
    for i_k=1:n_k
        plot(se_ratio,wc_prob_wald(plot_ellipse_M(i_M),:,i_k),linestyles{i_k},'LineWidth',1);
    end
    hold off;
    the_xlim = xlim;
    line(the_xlim,[1 1]*(1-signif),'Color','k','LineStyle','-'); % Mark nominal coverage level on vertical axis
    xlim(the_xlim);
    ylim([0 1]);
    label_axis(false);
    title(sprintf('%s%3.2g%s','$M=',M(plot_ellipse_M(i_M)),'$'),'Interpreter','Latex','FontSize',14);
    if i_M==1
        label_legend(k,'k');
    end
end


%% Plotting functions

function label_axis(univariate)
    if univariate
        xlabel('$\sqrt{\mathrm{aVar}(\hat{\delta}_h)/\mathrm{aVar}(\hat{\beta}_h)}$','Interpreter','Latex');
    else
        xlabel('$\sqrt{\lambda_{\mathrm{min}}(\mathrm{aVar}(\hat{\mbox{\boldmath$\delta$}})\mathrm{aVar}(\hat{\mbox{\boldmath$\beta$}})^{-1})}$','Interpreter','Latex');
    end
    set(gca,'FontSize',12);
end

function label_legend(val,s)
    legend(strcat(strcat(s,'='), arrayfun(@num2str, val, 'UniformOutput', 0)),'Location','SouthEast','FontSize',12);
end