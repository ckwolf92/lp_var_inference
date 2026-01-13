% setup_appl()   Setup application object.
% Input:
% - name (char):    Name of application.
% - df_raw (table): Imported Excel table from Ramey (2016).
% - est (struct):   Estimation settings
% - horzs (vector): IRF horizons.
% - p (numeric): Lag length.
% - fun_clean: Function handle for application-specific cleaning script.

function appl = setup_appl(name, df_raw, est, horzs, p, fun_clean)

    appl             = struct;
    appl.name        = name;
    appl.data.df_raw = df_raw;

    [appl.data.y, appl.data.shock, appl.data.date,...
        appl.data.yname, appl.data.shockname]     = fun_clean(df_raw);

    appl.est                = est;
    appl.est.resp_ind       = 2:length(appl.data.yname)+1;
    appl.est.horzs          = horzs;
    appl.est.p              = p;

end
