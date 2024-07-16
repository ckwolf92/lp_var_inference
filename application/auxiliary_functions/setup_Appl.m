% setup_Appl()   Setup application object.
% Input:
% - name (char):    Name of application.
% - df_raw (table): Imported Excel table from Ramey (2016).
% - est (struct):   Estimation settings
% - horzs (vector): IRF horizons.
% - p (numeric): Lag length.
% - fun_clean: Function handle for application-specific cleaning script.

function Appl = setup_Appl(name, df_raw, est, horzs, p, fun_clean)

    Appl             = struct;
    Appl.name        = name;
    Appl.data.df_raw = df_raw;

    [Appl.data.y, Appl.data.shock, Appl.data.date,...
        Appl.data.yname, Appl.data.shockname]     = fun_clean(df_raw);
    Appl.est                = est;
    Appl.est.resp_ind       = 2:length(Appl.data.yname)+1;
    Appl.est.horzs          = horzs;
    Appl.est.p              = p;


end
