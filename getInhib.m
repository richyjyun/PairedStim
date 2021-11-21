%% Get inhibitory response 
% First spike after stimulation that is not an evoked spike (greater than
% dt)

function inhib = getInhib(dt,spikes,stim)

mintime = max(dt);
temp = discretize(spikes,stim);
allinds = unique(temp); allinds(isnan(allinds)) = [];
inhib = nan(1,length(stim));
for i = 1:length(allinds)
    times = spikes(temp==allinds(i))-stim(allinds(i));
    ind = find(times>mintime,1);
    
    if(~isempty(ind))
        inhib(allinds(i)) = times(ind);
    end
end


end

