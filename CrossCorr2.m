function [cor,lags] = CrossCorr2(ts1, ts2, bin, cortime, plt)
%% Calculates cross correlation between two timestamps
% ts1 and ts2 are the timestamps
% fs is the sampling rate
% bin is the binsize for the crosscorrelation
% cortime is the maximum lag time

maxtime = max(max(ts1),max(ts2));
corvec_1 = zeros(round(maxtime/bin),1);
corvec_1(round(ts1/bin)) = 1;

corvec_2 = zeros(round(maxtime/bin),1);
corvec_2(round(ts2/bin)) = 1;

[cor,lags] = xcorr(corvec_2,corvec_1,round(cortime/bin));
if(all(corvec_1==corvec_2))
    cor(lags==0) = 0;
end

cor = cor/length(ts1);

lags = lags*bin;

if(exist('plt') && plt == 1)
    figure; bar(lags,cor,1,'k'); yl = ylim; ylim([0,yl(2)]);
    xlabel('Seconds'); ylabel('Probability');
end

end