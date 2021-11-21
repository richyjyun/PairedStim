%% Get evoked spikes from spike and stimulus timestamp data

function [ES, dt, prc, norm] = getESpikes(spikes, stim)

%% Setup
max_eSpike_time = 0.02;
min_eSpike_time = 0.0010;

bin = 0.0005; cortime = 0.2;
[cor,lags] = CrossCorr2(stim, spikes, bin, cortime);

norminds = find(lags>=-0.1,1):find(lags>=-0.01,1);
norm = cor(norminds);

zerolag = find(lags==0);

smoothwin = 5; % window to smooth PSTH
PSTH = smoothdata(cor,'gaussian',smoothwin); % smooth

baselim = -0.02;
stimfreq = length(stim)/(stim(end)-stim(1));
if(stimfreq > 30), baselim = -0.01; end
basewin = [find(lags==baselim), find(lags==-0.002)]; % window for defining baseline

%% Find peaks
% set threshold for determining evoked spikes. 2*std is flexible.
thresh = mean(cor(basewin)) + std(cor(basewin));

% find peaks
[p,locs] = findpeaks(PSTH(zerolag:(zerolag+round(max_eSpike_time/bin))),'MinPeakHeight',thresh);
p(lags(locs+zerolag-1)<min_eSpike_time) = [];
locs(lags(locs+zerolag-1)<min_eSpike_time) = [];

locs = locs(find(p==max(p)));

% if no peaks above threshold, there are no evoked spikes.
if(isempty(locs)), ES = nan; dt = nan; prc = nan; return; end

%% Find right limit
% find minimum in PSTH from the last PSTH peak to 30 ms.
leftlim = locs(end)+zerolag-1;
rightlim = find(lags>=(max(lags(locs+zerolag-1))+0.01),1);
[~,trough1] = min(PSTH(leftlim:rightlim));
trough1 = trough1 + leftlim - 1;

% find the first instance PSTH falls below threshold after the last espike peak.
% -2*std is flexible
minthresh = mean(cor(basewin)) - std(cor(basewin));
trough2 = find(PSTH(leftlim:rightlim)<minthresh,1);
trough2 = trough2 + leftlim;
if(isempty(trough2)), trough2 = rightlim; end

% the earliest of the previous two metrics is when espikes end
righttrough = min(trough1,trough2);
max_eSpike_time = lags(righttrough);

%% Find left limit
% find minimum in PSTH from 1 ms to the first PSTH peak
leftlim = find(lags>=0.001,1);
rightlim = locs(1)+zerolag-1;
[~,trough1] = min(PSTH(leftlim:rightlim));
trough1 = trough1 + leftlim - 1;

% find the last instance PSTH falls below threshold before the first espike peak.
% -2*std is flexible
trough2 = find(PSTH(leftlim:rightlim)<minthresh,1,'last');
trough2 = trough2 + leftlim - 2;

% the latest of the previous two metrics is when espikes end
lefttrough = max([trough1,trough2]);
min_eSpike_time = max(lags(lefttrough),min_eSpike_time);

if(isempty(max_eSpike_time) || isempty(min_eSpike_time))
    ES = nan; dt = nan; prc = nan;
    return;
end

%% Define evoked spikes
inds = discretize(spikes,stim);
inds(isnan(inds)) = 1;
spike_stim_diff = spikes - stim(inds);
ES = spike_stim_diff >= min_eSpike_time & spike_stim_diff <= max_eSpike_time;

%% Delay of each spike from stim time
dt = spike_stim_diff(ES);

%% Percentage of evoked spikes
prc = sum(ES)/length(stim)*100;





