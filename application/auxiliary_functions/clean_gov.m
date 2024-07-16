function [y, shock, date, yname, shockname] = clean_gov(data_raw)

df      = data_raw;
df.date = clean_date(df.quarter);

% -------------------------------------------------------------------------
% CONSTRUCT VARIABLES
% -------------------------------------------------------------------------

% define quartic trends
df.t  = (1:height(df))'; 
df.t2 = df.t.^2;
df.t3 = df.t.^3;
df.t4 = df.t.^4;

df.taxrate     = df.nfedreceipts./ df.ngdp;  % average tax rate 
df.ncndsv      = df.ncnd + df.ncsv;          % nominal nondurable plus services consumption
df.tothourspc  = df.tothours ./ df.pop;      % per capita hours worked
df.rwbus       = df.nwbus./ df.pbus;         % real compensation in business, divided by deflator for business
df.lrwbus      = log(df.rwbus);
df.ltothourspc = log(df.tothourspc);
df.lpgdp       = log(df.pgdp);
df.infl        = [NaN; 400*log(df.pgdp(2:end) ./ df.pgdp(1:end-1))]; % inflation
df.rint        = df.tbill3 - df.infl;  % real interest rate


% Define real NIPA variables with common GDP deflator
df.rgdp         = 100*df.ngdp         ./ df.pgdp;
df.rgov         = 100*df.ngov         ./ df.pgdp;
df.rfedreceipts = 100*df.nfedreceipts ./ df.pgdp;
df.rcons        = 100*df.ncons        ./ df.pgdp;
df.rtotinv      = 100*df.ntotinv      ./ df.pgdp;
df.rcndsv       = 100*df.ncndsv       ./ df.pgdp;
df.rcdur        = 100*df.ncdur        ./ df.pgdp;
df.rnri         = 100*df.nnri         ./ df.pgdp;
df.rres         = 100*df.nres         ./ df.pgdp;


df.lrgdp = log(df.rgdp);
df.lrgov = log(df.rgov);
df.lrtax = log(df.rfedreceipts);


% estimate trend GDP for Gordon-Krenn normalization;
controls         = [ones(height(df), 1), df.t, df.t2];
df.lyquad = controls*inv(controls'*controls)*controls' * df.lrgdp;
df.yquad  = exp(df.lyquad);


% Ramey narrative military news, divided by normalization 
df.newsy = [NaN; 100*df.rameynews(2:end) ./ (df.yquad(1:end-1) .* df.pgdp(1:end-1))]; 

% normalize a la Gordon-Krenn

df.rgdpx         = df.rgdp         ./ df.yquad;
df.rgovx         = df.rgov         ./ df.yquad;
df.rfedreceiptsx = df.rfedreceipts ./ df.yquad;
df.rcndsvx       = df.rcndsv       ./ df.yquad;
df.rcdurx        = df.rcdur        ./ df.yquad;
df.rnrix         = df.rnri         ./ df.yquad;
df.rresx         = df.rres         ./ df.yquad;
df.rconsx        = df.rcons        ./ df.yquad;
df.rtotinvx      = df.rtotinv      ./ df.yquad;
 
 
 
% RENAME TO SHORTER NAMES;  
df = renamevars(df, {'rgdpx'   , 'rgovx'      , 'rconsx', 'rtotinvx', ...
                 'taxrate', 'rcndsvx'    , 'rcdurx', 'rnrix', ...
                 'rresx'  , 'ltothourspc', 'lrwbus'}, ...
                {'y' ,'g'  ,'c'   ,'inv', ...
                'tax','cns','cdur','nri', ...
                'res','hr' ,'w'});

% 
keep     = find(year(df.date) == 1947 & month(df.date) ==6):find(year(df.date) == 2013 & month(df.date) ==12);
df       = df(keep, :);
controls = [ones(height(df), 1), df.t, df.t2];


M     = eye(size(controls, 1)) - controls*inv(controls'*controls)*controls';
y  = M*[df.y, df.g, df.tax];
shock = M*df.newsy;

% save
date      = df.date;
yname      = {'y', 'g', 'tax'};
shockname = 'newsy';

end