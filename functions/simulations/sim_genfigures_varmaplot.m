% sim_genfigures_varmaplot.m  VARMA: Select DGPs, CI procedures, and line specs for
%                             sim_genfigures.m

switch dgp_type(7:8)
    case 'fi'
        p_sel = [1,2,4,8];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'es'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'wo'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'LP'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            2 1];
        line_colors = [204/255 0/255 0/255; 102/255 178/255 255/255];
        line_specs = {'-', '-.'};
    case 'bf'
        p_sel = [1,2,4,8];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            1 4;
            2 1;
            2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
    case 'be'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            1 4;
            2 1;
            2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
    case 'bw'
        p_sel = [1];
        rhos_sel   = [0.9];
        proc_names = {'VAR', 'VAR$_b$', 'LP', 'LP$_b$'};
        procs = [1 1; % first index: inference procedure; second index: type of confidence interval
            1 4;
            2 1;
            2 4];
        line_colors = [204/255 0/255 0/255; 0.5 * [1 1 1] + 0.5 * [204/255 0/255 0/255]; ...
            102/255 178/255 255/255; 0.5 * [1 1 1] + 0.5 * [102/255 178/255 255/255]];
        line_specs = {'-', '-.', ':', '--'};
end