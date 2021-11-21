%% Analyzing coherence between the Pre and Post sites of paired stimulation

%% Load data
if(~exist('SL'))
    SLfile = 'C:\SL.mat';
    load(SLfile);
end


%% Extract all variables
PreCoh = {}; PostCoh = {}; Cohf = {};
% Need to loop, can't use extractfield
for i = 1:length(SL)
    PreCoh{i} = SL(i).PreCoh; 
    PostCoh{i} = SL(i).PostCoh;
    Cohf{i} = SL(i).Cohf;
end
temp = cellfun(@isempty,Cohf);
temp = find(~temp,1);
Cohf = Cohf{temp};

%% Average differences in coherence
temp = cellfun(@(x,y) smooth(x,20)-smooth(y,20), PostCoh, PreCoh, 'Uniformoutput', false);
% Convert to matrix. have to loop because of empty ones
N = max(cellfun(@length,PreCoh));
Cohdiff = nan(N,length(PreCoh));
for i = 1:length(PreCoh)
    if(~isempty(temp{i}))
        Cohdiff(:,i) = temp{i};
    end
end

% Averaged across groups. Run cs_AnalyzeEvokedSpikes before this
figure;
for i = 1:4
    hold on; plot(Cohf,nanmean(Cohdiff(:,ind{i}),2),colors{i});
end
xlim([0,100])

% Within frequency bands (all)
bands = {'theta','alpha','beta','low gamma','mid gamma','high gamma'};
freqs = [4,8; 8,12; 15,30; 30,50; 50,80; 80,150];
figure; pvalCoh = [];
for i = 1:length(bands)
    low = find(Cohf>=freqs(i,1),1);
    high = find(Cohf>=freqs(i,2),1);
    temp = Cohdiff(low:high,:);
    hold on; boxplot(nanmean(temp,2), 'position', i, 'notch', 'on', 'symbol', 'w')
    [~,pvalCoh(i)] = ttest(nanmean(temp,2));
end

% Frequency bands within groups
% Within frequency bands (all)
bands = {'theta','alpha','beta','low gamma','mid gamma','high gamma'};
freqs = [4,7; 8,12; 15,30; 30,50; 50,80; 80,200];
figure; 
for i = 1:4
    subplot(1,4,i)
    difftemp = Cohdiff(:,ind{i});
    p = [];
    for j = 1:length(bands)
        low = find(Cohf>=freqs(j,1),1);
        high = find(Cohf>=freqs(j,2),1);
        temp = difftemp(low:high,:);
        hold on; boxplot(nanmean(temp), ...
            'position', j, 'notch', 'on', 'symbol', 'w', 'width', 0.5)
        [p(j),~] = signrank(nanmean(temp));
    end
    sig = find(p<0.05);
    sigstar(mat2cell(sig,1,ones(1,length(sig))),p(sig),0,20,0);
    ylim([-0.1,0.05])
end


figure;
difftemp = Cohdiff(:,ind{i});
p = [];
for j = 1:length(bands)
    subplot(1,length(bands),j)
    
    low = find(Cohf>=freqs(j,1),1);
    high = find(Cohf>=freqs(j,2),1);
    temp = difftemp(low:high,:);
    plot(Cohf(low:high),temp);
end







