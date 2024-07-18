function [y, shock, date, yname, shockname] = clean_tech(data_raw)

df       = data_raw;

% read and adjust data
df.date = clean_date(df.quarter);
keep    = find(year(df.date)==1949 & month(df.date)==6):...
         find(year(df.date) == 2009 & month(df.date) == 12);
df = df(keep, :);

df.t = (1:height(df))';
df.t2 = df.t.^2;

df.xtot    = df.rgdp     ./ df.tothours;
df.rstockp = df.stockp_sh./ df.pgdp;

for var = {'rgdp' 'rcons' 'rnri' 'ybus' 'tothours' 'nbus' 'rstockp'}
    df{:, var} = 1000*df{:, var} ./ df.pop;
end

for var = {'rgdp' 'rcons' 'rnri' 'xtot' 'ybus' 'xbus' 'nbus' 'tothours' 'pgdp' 'rstockp'}
    df{:, strcat('l', var)} = log(df{:, var});
    df{:, strcat('dl', var)} = [NaN; df{2:end, strcat('l', var)} - df{1:end-1, strcat('l', var)}];
end

for var = {'ltfp' 'ltfp_util' 'ltfp_I' 'ltfp_I_util'}
    df{:, strcat('d', var)} = 400*[NaN; df{2:end, var} - df{1:end-1,  var}];
end

controls = [ones(height(df),1), df.t, df.t2];
M        = eye(size(controls, 1)) - controls*inv(controls'*controls)*controls';
y        = M*[df.lrgdp, df.lrstockp, df.lxtot];
shock    = M* df.ford_tfp; 

% save
date      = df.date;
yname      = {'lrgdp', 'lrstockp', 'lxtot'};
shockname = 'ford_tfp';

end
