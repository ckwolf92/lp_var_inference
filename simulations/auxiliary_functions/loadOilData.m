% load baseline data

load OilDataM
sampleDates    = datesStringRaw(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)));

% load proxy

ncontract = 15;       % select futures contract 

load OilSurprisesMLog

proxyRaw = [oilProxiesWTIM(:,ncontract)]; 

smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));

smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));

proxy = proxyRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR

[T,np] = size(proxy);
k = 1; % index of variable(s) to be instrumented

