function [Dist,windowRange,methodUsed] = localPeakValleyDistinguishability_v2(lineProfile,Nratio,method,dT,windowSec,peakSpreadThreshold)
if nargin<3 || isempty(method), method = 'original'; end
if nargin<4, dT = []; end
if nargin<5 || isempty(windowSec), windowSec = 0.1; end
if nargin<6 || isempty(peakSpreadThreshold), peakSpreadThreshold = 0.1; end

x = lineProfile(:);
x = x(isfinite(x));
method = lower(strrep(method,'-',''));
windowRange = [1 numel(x)];

if ismember(method,{'stablewindow','slidingaverage','windowaverage'})
    Dist = meanSlidingPVD(x,Nratio,dT,windowSec);
    methodUsed = "sliding";
    return;
end

if ismember(method,{'adaptive','auto'})
    [~,~,peakVals] = peakValleyStats(x,Nratio);
    peakSpreadRatio = spreadRatio(peakVals);
    if peakSpreadRatio>peakSpreadThreshold
        Dist = meanSlidingPVD(x,Nratio,dT,windowSec);
        methodUsed = "sliding";
    else
        [deltaVal,pooledStd] = peakValleyStats(x,Nratio);
        Dist = safeRatio(deltaVal,pooledStd);
        methodUsed = "original";
    end
    return;
end

[deltaVal,pooledStd] = peakValleyStats(x,Nratio);
switch method
    case {'original','classic','tailstd'}
        denomVal = pooledStd;
        methodUsed = "original";
    case {'modified','linenormalized'}
        denomVal = pooledStd/std(x,0,'omitnan');
        methodUsed = "modified";
    otherwise
        error('Unknown PVD method: %s',method);
end

Dist = safeRatio(deltaVal,denomVal);
end

function Dist = meanSlidingPVD(x,Nratio,dT,windowSec)
lineN = numel(x);
if lineN<10
    Dist = NaN;
    return;
end

[windowN,stepN] = windowSamples(lineN,dT,windowSec);
starts = 1:stepN:(lineN-windowN+1);
if starts(end)~=(lineN-windowN+1)
    starts = [starts lineN-windowN+1];
end

windowPVD = nan(numel(starts),1);
for ii = 1:numel(starts)
    currX = x(starts(ii):(starts(ii)+windowN-1));
    [deltaVal,pooledStd] = peakValleyStats(currX,Nratio);
    windowPVD(ii) = safeRatio(deltaVal,pooledStd);
end

Dist = mean(windowPVD,'omitnan');
end

function [deltaVal,pooledStd,peakVals,valleyVals] = peakValleyStats(x,Nratio)
sampleN = max(1,round(numel(x)*Nratio));
peakVals = maxk(x,sampleN);
valleyVals = mink(x,sampleN);
deltaVal = mean(peakVals)-mean(valleyVals);
pooledStd = sqrt((var(peakVals,0)+var(valleyVals,0))/2);
end

function val = spreadRatio(vals)
meanVal = abs(mean(vals));
stdVal = std(vals,0);
if meanVal==0
    val = Inf;
else
    val = stdVal/meanVal;
end
end

function val = safeRatio(numVal,denomVal)
if denomVal==0 || isnan(denomVal)
    val = NaN;
else
    val = numVal/denomVal;
end
end

function [windowN,stepN] = windowSamples(lineN,dT,windowSec)
if ~isempty(dT) && isfinite(dT(1)) && dT(1)>0
    windowN = max(3,round(windowSec/dT(1)));
    windowN = min(windowN,lineN);
    stepN = max(1,round(windowN/5));
else
    windowN = min(max(10,round(0.45*lineN)),lineN);
    stepN = max(1,round(0.05*lineN));
end
end
