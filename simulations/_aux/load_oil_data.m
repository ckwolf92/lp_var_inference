% load baseline data

load OilDataM

data = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100];
data = data(169:end,:);

% load proxy

load OilSurprisesMLog

proxy = oilProxiesWTIM(:,15);

