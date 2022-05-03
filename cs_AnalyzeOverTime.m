%% Run cs_AnalyzeEvokedSpikes and cs_AnalyzeCCEPs first

%% Everything grouped
figure; comparisontype = 'LSD'; yl = [-3.5,3.5];
maxcompare = 1; keep = Cond;

% EVOKED SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRE
PreOverTime = {};
temp = min([max(cellfun(@length, PreOverTime)),60]);
bw = 10; x = linspace(bw/2,temp*bw-bw/2,temp);
for i = 1:length(PreES)
    if(isnan(PreES{i})), continue; end
    bins = 0:bw:(PreStim{i}(end));
    Stimcount = histcounts(PreStim{i},bins);
    EScount = histcounts(PreSpikes{i}(PreES{i}),bins);
    PreOverTime{i} = EScount./Stimcount;
end

% Convert to matrix
PreAll = nan(temp,length(PreOverTime));
for i = 1:length(PreOverTime)
    if(isempty(PreOverTime{i})),continue; end
    inds = 1:min(length(PreOverTime{i}),temp);
    PreAll(inds,i) = PreOverTime{i}(inds);
end

% zscore
PreAllES = (PreAll - nanmean(PreAll)) ./ nanstd(PreAll);

% Plot
subplot(3,2,1);
PlotOverTime(x,PreAllES(:,keep),comparisontype,'Evoked spike prob',yl,1,maxcompare);
title('Pre-test');

% POST
PostOverTime = {};
bw = 10; 
for i = 1:length(PostES)
    if(isnan(PostES{i})), continue; end
    bins = 0:bw:(PostStim{i}(end));
    Stimcount = histcounts(PostStim{i},bins);
    EScount = histcounts(PostSpikes{i}(PostES{i}),bins);
    PostOverTime{i} = EScount./Stimcount;
end

% Convert to matrix
temp = min([max(cellfun(@length, PostOverTime)),60]);
PostAll = nan(temp,length(PostOverTime));
for i = 1:length(PostOverTime)
    if(isempty(PostOverTime{i})),continue; end
    inds = 1:min(length(PostOverTime{i}),temp);
    PostAll(inds,i) = PostOverTime{i}(inds);
end
PostOverTimeES = PostAll;

% zscore
PostAllES = (PostAll - nanmean(PostAll)) ./ nanstd(PostAll);

% Plot
subplot(3,2,2);
PlotOverTime(x,PostAllES(:,keep),comparisontype,'',yl,0,maxcompare);
title('Post-test');

% INHIBITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRE
PreOverTime = {};
bw = 10; 
for i = 1:length(PreInhib)
    bins = 0:bw:(PreStim{i}(end));
    Stimind = discretize(PreStim{i},bins);
    PreOverTime{i} = arrayfun(@(x) nanmedian(PreInhib{i}(Stimind==x)), 1:max(Stimind));
end

% Convert to matrix
temp = min([max(cellfun(@length, PreOverTime)),60]);
PreAll = nan(temp,length(PreOverTime));
for i = 1:length(PreOverTime)
    if(isempty(PreOverTime{i})),continue; end
    inds = 1:min(length(PreOverTime{i}),temp);
    PreAll(inds,i) = PreOverTime{i}(inds);
end

% zscore
PreAllIH = (PreAll - nanmean(PreAll)) ./ nanstd(PreAll);
 
% Plot
subplot(3,2,3);
PlotOverTime(x,PreAllIH(:,keep),comparisontype,'Ihibition duration',yl,0,maxcompare);

% POST
PostOverTimeIH = {};
bw = 10;
for i = 1:length(PostInhib)
    bins = 0:bw:(PostStim{i}(end));
    Stimind = discretize(PostStim{i},bins);
    PostOverTimeIH{i} = arrayfun(@(x) nanmedian(PostInhib{i}(Stimind==x)), 1:max(Stimind));
end

% Convert to matrix
temp = min([max(cellfun(@length, PostOverTimeIH)),60]);
PostAll = nan(temp,length(PostOverTimeIH));
for i = 1:length(PostOverTimeIH)
    if(isempty(PostOverTimeIH{i})),continue; end
    inds = 1:min(length(PostOverTimeIH{i}),temp);
    PostAll(inds,i) = PostOverTimeIH{i}(inds);
end

PostOverTimeIH = PostAll;

% zscore
PostAllIH = (PostAll - nanmean(PostAll)) ./ nanstd(PostAll);

% Plot
subplot(3,2,4);
PlotOverTime(x,PostAllIH(:,keep),comparisontype,'',yl,0,maxcompare);


% CCEPs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRE
PreEPTimep2p = {}; PostEPTimep2p = {}; 
for i = 1:length(PreEP)
    [trough,peak,t_ind,p_ind] = getEPMag(Range{i}/fs,PreEPTime{i},fs);
    PreEPTimep2p{i} = peak-trough; 
    
    [trough,peak,t_ind,p_ind] = getEPMag(Range{i}/fs,PostEPTime{i},fs);
    PostEPTimep2p{i} = peak-trough; 
end

% Consolidate Pre into matrix
N = max(cellfun(@length, PreEPTimep2p));
PreMag = nan(length(PreEPTimep2p),N);
for i = 1:length(PreEPTimep2p)
    PreMag(i,1:length(PreEPTimep2p{i})) = PreEPTimep2p{i};
end

% zscore
PreMagZ = (PreMag - nanmean(PreMag,2))./nanstd(PreMag,[],2);
PreMagZ = PreMagZ';

% Plot
subplot(3,2,5);
PlotOverTime(x,PreMagZ(:,keep),comparisontype,'CCEP magnitude',yl,0,maxcompare);

% POST
N = max(cellfun(@length, PostEPTimep2p));
PostMag = nan(length(PostEPTimep2p),N);
for i = 1:length(PostEPTimep2p)
    PostMag(i,1:length(PostEPTimep2p{i})) = PostEPTimep2p{i};
end

% zscore
PostMagZ = (PostMag - nanmean(PostMag,2))./nanstd(PostMag,[],2);
nancount = isnan(PostMagZ); 
nancount = [ones(size(nancount,1),1),nancount];
nancount = diff(nancount,[],2); 
shiftcount = [];
for i = 1:size(nancount,1)
    shiftcount(i) = find(nancount(i,:)==-1,1)-1;
    PostMagZ(i,:) = circshift(PostMagZ(i,:),-shiftcount(i));
end
PostMagZ = PostMagZ(:,1:60)';

% Plot
subplot(3,2,6);
PlotOverTime(x,PostMagZ(:,keep),comparisontype,'',yl,0,maxcompare);


%% By regions
% Separate the regions
ind = {};
change = ESchange; 
keep = Cond;
ind{1} = normISI < 0 & normISI > -4 & keep & change < 0;%normISI < 0 & ~Bad;
ind{2} = normISI > 0 & normISI<4 & change > 0 & keep;%normISI < 1 & normISI > 0 & ESchange > 0 & ~Bad;
ind{3} = normISI > 0 & normISI<4 & change < 0 & keep;%normISI > 0 & normISI < 10 & ESchange < 0 & ~Bad;
ind{4} = ~ind{1} & ~ind{2} & ~ind{3} & keep;

temp = min([max(cellfun(@length, PreOverTime)),60]);
x = linspace(bw/2,temp*bw-bw/2,temp);
figure; leg = 0; ylab = '';
for i = 1:length(ind)
    subplot(3,length(ind),i);
    
    temp = PostAllES(:,ind{i});
    if(i==1)
        leg=1; ylab= 'Evoked spike prob';
    else
        leg=0; ylab = '';
    end
    PlotOverTime(x,temp,comparisontype,ylab,yl,leg,maxcompare);
    
end

for i = 1:length(ind)
    subplot(3,length(ind),i+length(ind));
    
    temp = PostAllIH(:,ind{i});
    if(i==1)
        ylab= 'Inhibition duration';
    else
        ylab = '';
    end
    PlotOverTime(x,temp,comparisontype,ylab,yl,leg,maxcompare);
    
end

for i = 1:length(ind)
    subplot(3,length(ind),i+2*length(ind));
    
    temp = PostMagZ(:,ind{i});
    if(i==1)
        ylab= 'CCEP magnitude';
    else
        ylab = '';
    end
    PlotOverTime(x,temp,comparisontype,ylab,yl,leg,maxcompare);
    
end

%% Scatter plots
% Evoked spikes
FirstES = (nanmean(PostOverTimeES(1:3,:))-PreESProb)./PreESProb*100;
LastES = (nanmean(PostOverTimeES(end-2:end,:))-PreESProb)./PreESProb*100;
fixind = isnan(LastES);
LastES(fixind) = (nanmean(PostOverTimeES(28:30,fixind))-PreESProb(fixind))./PreESProb(fixind)*100;
ESdiff = LastES-FirstES;

% Raw ISI
figure;
for i = 1:3
    subplot(3,3,i);
    switch i
        case 1
            data = FirstES;
        case 2
            data = LastES;
        otherwise
            data = ESdiff;
    end
    colors = get(gca,'colororder');
    scatter(normISI(Cond & monkeyJ),data(Cond & monkeyJ),15,'o',...
        'MarkerEdgeColor',colors(1,:),'linewidth',1.5); hold on;
    scatter(normISI(Cond & monkeyL),data(Cond & monkeyL),15,'^',...
        'MarkerEdgeColor',colors(2,:),'linewidth',1.5);
    ylim([-100,100]); xlim([-30,30]); xl = xlim; yl = ylim;
    hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');
    xlabel('ISI (a.u.)'); ylabel('% Change');
    set(gca,'FontSize',10);
end

% Inhibition
PreIH = cellfun(@nanmedian, PreInhib);
FirstIH = (nanmean(PostOverTimeIH(1:3,:))-PreIH)./PreIH*100;
LastIH = (nanmean(PostOverTimeIH(end-2:end,:))-PreIH)./PreIH*100;
fixind = isnan(LastIH);
LastIH(fixind) = (nanmean(PostOverTimeIH(28:30,fixind))-PreIH(fixind))./PreIH(fixind)*100;
IHdiff = LastIH-FirstIH;

% Raw ISI
for i = 1:3
    subplot(3,3,i+3);
    switch i
        case 1
            data = FirstIH;
        case 2
            data = LastIH;
        otherwise
            data = IHdiff;
    end
    colors = get(gca,'colororder');
    scatter(normISI(Cond & monkeyJ),data(Cond & monkeyJ),15,'o',...
        'MarkerEdgeColor',colors(1,:),'linewidth',1.5); hold on;
    scatter(normISI(Cond & monkeyL),data(Cond & monkeyL),15,'^',...
        'MarkerEdgeColor',colors(2,:),'linewidth',1.5);
    ylim([-100,100]); xlim([-30,30]); xl = xlim; yl = ylim;
    hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');
    xlabel('ISI (a.u.)'); ylabel('% Change');
    set(gca,'FontSize',10);
end

% CCEP
if ~(size(PostMag,1) < 100)
    nancount = isnan(PostMag);
    nancount = [ones(size(nancount,1),1),nancount];
    nancount = diff(nancount,[],2);
    shiftcount = [];
    for i = 1:size(nancount,1)
        shiftcount(i) = find(nancount(i,:)==-1,1)-1;
        PostMag(i,:) = circshift(PostMag(i,:),-shiftcount(i));
    end
    PostMag = PostMag(:,1:60)';
end

FirstEP = (nanmean(PostMag(1:3,:))-PreEPp2p)./PreEPp2p*100;
LastEP = (nanmean(PostMag(end-2:end,:))-PreEPp2p)./PreEPp2p*100;

firstnan = find(isnan(FirstEP));
lastnan = find(isnan(LastEP));
for i = 1:length(lastnan)
    if(any(lastnan(i)==firstnan)), continue; end
        temp = PostMag(:,lastnan(i));
        goodind = find(~isnan(temp),1,'last');
        LastEP(lastnan(i)) = (nanmean(temp(goodind-2:goodind))-PreEPp2p(lastnan(i)))...
            ./PreEPp2p(lastnan(i))*100;
end

EPdiff = LastEP-FirstEP;

% Raw ISI
for i = 1:3
    subplot(3,3,i+6);
    switch i
        case 1
            data = FirstEP;
        case 2
            data = LastEP;
        otherwise
            data = EPdiff;
    end
    colors = get(gca,'colororder');
    scatter(normISI(Cond & monkeyJ),data(Cond & monkeyJ),15,'o',...
        'MarkerEdgeColor',colors(1,:),'linewidth',1.5); hold on;
    scatter(normISI(Cond & monkeyL),data(Cond & monkeyL),15,'^',...
        'MarkerEdgeColor',colors(2,:),'linewidth',1.5);
    ylim([-100,100]); xlim([-30,30]); xl = xlim; yl = ylim;
    hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');
    xlabel('ISI (a.u.)'); ylabel('% Change');
    set(gca,'FontSize',10);
end


%% Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












