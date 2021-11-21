%% Analyze stimulus responses including evoked spikes and inhibition

%% Load data
if(~exist('SL'))
    SLfile = 'C:\SL.mat';
    load(SLfile);
end

%% Extract all variables
ISI = extractfield(SL,'ISI');
Animal = extractfield(SL,'Animal');
PreStim = {}; PostStim = {};
PreSpikes = {}; PostSpikes = {};
PreES = {}; PostES = {}; 
Predt = {}; Postdt = {};
PreInhib = {}; PostInhib = {};
PreNorm = {}; PostNorm = {};
for i = 1:length(SL)
    
    SD = SL(i);
    
    % Set variables
    PreStim{i} = SD.PreStim; PostStim{i} = SD.PostStim;
    PreSpikes{i} = SD.PreSpikes; PostSpikes{i} = SD.PostSpikes;
    
    % Extract evoked spikes
    [PreES{i}, Predt{i}, ~, PreNorm{i}] = getESpikes(SD.PreSpikes, SD.PreStim);
    [PostES{i}, Postdt{i}, ~, PostNorm{i}] = getESpikes(SD.PostSpikes, SD.PostStim);
    
    % Extract inhibition duration
    PreInhib{i} = getInhib(Predt{i},SD.PreSpikes, SD.PreStim);    
    PostInhib{i} = getInhib(Postdt{i},SD.PostSpikes, SD.PostStim);
    
end

% Set indices for analysis
Bad = cell2mat(extractfield(SL,'Bad'));
Control = extractfield(SL,'Control');
Cond = ~(Bad|~isnan(Control));

% Set indices for animals
monkeyJ = cellfun(@(x) strcmp(x,'Jafar'), Animal);
monkeyL = cellfun(@(x) strcmp(x,'Lorde'), Animal);

%% Evoked spikes
%% Statistics on evoked spike probabilities - Fisher's exact test
skip = find(Bad);
pvalES = nan(1,length(PreES));
for i = 1:length(PreES)
    
    if(any(intersect(skip,i))), continue; end
    
    n1 = sum(PreES{i}); N1 = length(PreStim{i});
    n2 = sum(PostES{i}); N2 = length(PostStim{i});
    
    if(isnan(n1) || isnan(n2)), continue; end
    
    [~,pvalES(i),~] = fishertest([n1,N1;n2,N2]);
    
end

%% Evoked spike probabilty changes vs ISI
% Calculate percent change of probabilities
PreESProb = cellfun(@(x,y) sum(x)/length(y), PreES, PreStim);
PostESProb = cellfun(@(x,y) sum(x)/length(y), PostES, PostStim);
ESchange = (PostESProb-PreESProb)./PreESProb * 100;

pvallim = 0.05;
sig = pvalES < pvallim;

% Raw ISI
figure; subplot(1,2,1); scatter(ISI(Cond & monkeyJ & ~sig),ESchange(Cond & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(ISI(Cond & monkeyL & ~sig),ESchange(Cond & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(ISI(Cond & monkeyJ & sig),ESchange(Cond & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(ISI(Cond & monkeyL & sig),ESchange(Cond & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k'); 
xlabel('ISI (ms)'); ylabel('% Change');
set(gca,'FontSize',12); title('Raw ISI')

legend({'Monkey J','Monkey L'});

% Normalized ISI
% Calculate normalized ISI
normISI = cellfun(@(x,y) x/1000/nanmedian(y), mat2cell(ISI,1,ones(1,length(ISI))), Predt);

subplot(1,2,2); scatter(normISI(Cond & monkeyJ & ~sig),ESchange(Cond & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(normISI(Cond & monkeyL & ~sig),ESchange(Cond & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(normISI(Cond & monkeyJ & sig),ESchange(Cond & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(normISI(Cond & monkeyL & sig),ESchange(Cond & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Normalized ISI')

legend({'Monkey J','Monkey L'});

% % Normalized ISIs and normalized ES
% normPreESProb = cellfun(@(x,y,z) sum(x)/length(y)/median(z), PreES, PreStim, PreNorm);
% normPostESProb = cellfun(@(x,y,z) sum(x)/length(y)/median(z), PostES, PostStim, PostNorm);
% normESchange = (normPostESProb-normPreES)./normPreESProb * 100;
% 
% subplot(1,3,3); scatter(normISI(~Bad & monkeyJ & ~sig),normESchange(~Bad & monkeyJ & ~sig),45,'k^','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyL & ~sig),normESchange(~Bad & monkeyL & ~sig),40,'ko','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyJ & sig),normESchange(~Bad & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyL & sig),normESchange(~Bad & monkeyL & sig),45,'filled','ko','linewidth',1.5);
% ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
% hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
% xlabel('Normalized ISI'); ylabel('% Change');
% set(gca,'FontSize',12); title('Normalized ISI, Normalized ES')
% 
% legend({'Monkey J','Monkey L'});

%% Raw evoked spike change vs normalized changes
% 2 sample vs paired
% significantly different, but around 5% difference on average
figure; scatter(ESchange(Cond & monkeyJ),normESchange(Cond & monkeyJ),45,'k^');
hold on; scatter(ESchange(Cond & monkeyL),normESchange(Cond & monkeyL),40,'ko');
xlim([-100,100]); ylim([-100,100]);
xl = xlim; yl = xlim; hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');
mdl = fitlm(ESchange(Cond),normESchange(Cond));
b = mdl.Coefficients{1,1}; m = mdl.Coefficients{2,1};
hold on; plot(xl, m*xl+b, 'r--', 'linewidth', 1.5);
xlabel('Raw ES % Change'); ylabel('Normalized ES % Change');
[h,p] = ttest(ESchange(Cond),normESchange(Cond));
[p,h] = signrank(ESchange(Cond),normESchange(Cond));

%% Firing rate before and after conditioning
% Not significant
PreFR = cellfun(@(x,y) (length(x)-sum(y))/(x(end)-x(1)), PreSpikes, PreES);
PostFR = cellfun(@(x,y) (length(x)-sum(y))/(x(end)-x(1)), PostSpikes, PostES)+0.5;
figure; scatter(PreFR(Cond & monkeyJ),PostFR(Cond & monkeyJ),45,'k^');
hold on; scatter(PreFR(Cond & monkeyL),PostFR(Cond & monkeyL),40,'ko');
xlim([0,30]); ylim([0,30]); hold on; plot(xl,yl,'k');
xlabel('Pre FR (Hz)'); ylabel('Post FR (Hz)');

[h,p] = ttest2(PreFR(Cond),PostFR(Cond));
[p,h] = ranksum(ESchange(Cond),normESchange(Cond));

%% Initial evoked probability vs changes
% Not significant
PreESProb = cellfun(@(x,y) sum(x)/length(y), PreES, PreStim)*100;
figure; scatter(PreESProb(Cond & monkeyJ),ESchange(Cond & monkeyJ),45,'k^');
hold on; scatter(PreESProb(Cond & monkeyL),ESchange(Cond & monkeyL),40,'ko');
xlim([0,100]); ylim([-100,100]);
xl = xlim; yl = xlim; hold on; plot(xl,[0,0],'k');
mdl = fitlm(PreESProb(Cond),ESchange(Cond));
b = mdl.Coefficients{1,1}; m = mdl.Coefficients{2,1};
hold on; plot(xl, m*xl+b, 'r--', 'linewidth', 1.5);
xlabel('Basline ES Prob'); ylabel('Normalized ES % Change');

%% Inhibition
%% Statistics on inhibition duration - Wilcoxon ranksum
pvalIH = nan(1,length(PreInhib));
for i = 1:length(PreInhib)
       
    a = PreInhib{i}; a(isnan(a)) = [];
    b = PostInhib{i}; b(isnan(b)) = [];
   
    if(isempty(a) || isempty(b)), continue; end
    
    [pvalIH(i),~] = ranksum(a,b);
end

%% Inhibition duration change vs ISI
pvallim = 0.05;
sig = pvalIH < pvallim;

% Raw ISI with raw IH change
IHchange = cellfun(@(x,y)(nanmedian(x)-nanmedian(y))/nanmedian(x)*100, PostInhib, PreInhib);
figure; subplot(1,3,1); scatter(ISI(~Bad & monkeyJ & ~sig),IHchange(~Bad & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(ISI(~Bad & monkeyL & ~sig),IHchange(~Bad & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(ISI(~Bad & monkeyJ & sig),IHchange(~Bad & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(ISI(~Bad & monkeyL & sig),IHchange(~Bad & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-40,40]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k'); 
xlabel('Stim Time (ms)'); ylabel('% Change');
set(gca,'FontSize',12); title('Raw ISI, Inhib')

legend({'Monkey J','Monkey L'});

% Normalized ISI with raw IH change
subplot(1,2,2); scatter(normISI(~Bad & monkeyJ & ~sig),IHchange(~Bad & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL & ~sig),IHchange(~Bad & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyJ & sig),IHchange(~Bad & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL & sig),IHchange(~Bad & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-40,40]); xlim([-10,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized Stim Time'); ylabel('% Change');
set(gca,'FontSize',12); title('Normalized ISI, Inhib')

legend({'Monkey J','Monkey L'});

% % Normalized ISI with normalized IH change (not necessary, inhibition not
% % usually correlated with firing rate. 
% normIHchange = cellfun(@(x,y,a,b)(nanmedian(x)/nanmedian(b)-nanmedian(y)/nanmedian(a))/(nanmedian(x)/nanmedian(a))*100, PostInhib, PreInhib, PreNorm, PostNorm);
% subplot(1,3,3); scatter(normISI(~Bad & monkeyJ & ~sig),normIHchange(~Bad & monkeyJ & ~sig),45,'k^','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyL & ~sig),normIHchange(~Bad & monkeyL & ~sig),40,'ko','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyJ & sig),normIHchange(~Bad & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
% hold on; scatter(normISI(~Bad & monkeyL & sig),normIHchange(~Bad & monkeyL & sig),45,'filled','ko','linewidth',1.5);
% ylim([-40,40]); xlim([-10,25]); xl = xlim; yl = ylim; 
% hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
% xlabel('Normalized Stim Time'); ylabel('% Change');
% set(gca,'FontSize',12); title('Normalized ISI, Raw Inhib')
% 
% legend({'Monkey J','Monkey L'});

%% Evoked spike vs Inhibition
%% Evoked spike probability changes vs Inhibition duration changes
% Everything
figure; subplot(1,2,1);
scatter(ESchange(Cond & monkeyJ),IHchange(Cond & monkeyJ),45,'k^');
hold on; scatter(ESchange(Cond & monkeyL),IHchange(Cond & monkeyL),40,'ko');
ylim([-100,100]); xlim([-100,100])
xl = xlim; yl = ylim; hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');

mdl1 = fitlm(ESchange(Cond & monkeyJ),IHchange(Cond & monkeyJ)); 
b = mdl1.Coefficients{1,1}; m = mdl1.Coefficients{2,1};
hold on; plot(ESchange(Cond), m*ESchange(Cond)+b, 'r--', 'linewidth', 1.5);

mdl2 = fitlm(ESchange(Cond & monkeyL),IHchange(Cond & monkeyL)); 
b = mdl2.Coefficients{1,1}; m = mdl2.Coefficients{2,1};
hold on; plot(ESchange(Cond), m*ESchange(Cond)+b, 'b--', 'linewidth', 1.5);

mdl3 = fitlm(ESchange(Cond),IHchange(Cond));
b = mdl3.Coefficients{1,1}; m = mdl3.Coefficients{2,1};
hold on; plot(ESchange(Cond), m*ESchange(Cond)+b, 'k--', 'linewidth', 1.5);

% Normalized ES change. Statistically significant within each animal and
% while combined
subplot(1,2,2);
scatter(normESchange(Cond & monkeyJ),IHchange(Cond & monkeyJ),45,'k^');
hold on; scatter(normESchange(Cond & monkeyL),IHchange(Cond & monkeyL),40,'ko');
ylim([-100,100]); xlim([-100,100])
xl = xlim; yl = ylim; hold on; plot([0,0],yl,'k'); plot(xl,[0,0],'k');

mdl1 = fitlm(normESchange(Cond & monkeyJ),IHchange(Cond & monkeyJ)); 
b = mdl1.Coefficients{1,1}; m = mdl1.Coefficients{2,1};
hold on; plot(normESchange(Cond), m*normESchange(Cond)+b, 'r--', 'linewidth', 1.5);

mdl2 = fitlm(normESchange(Cond & monkeyL),IHchange(Cond & monkeyL)); 
b = mdl2.Coefficients{1,1}; m = mdl2.Coefficients{2,1};
hold on; plot(normESchange(Cond), m*normESchange(Cond)+b, 'b--', 'linewidth', 1.5);

mdl3 = fitlm(normESchange(Cond),IHchange(Cond));
b = mdl3.Coefficients{1,1}; m = mdl3.Coefficients{2,1};
hold on; plot(normESchange(Cond), m*normESchange(Cond)+b, 'k--', 'linewidth', 1.5);

%% Evoked spikes vs Inhibition, but split into sections of normalized ISI
% 2 stories:
% 1. Section 1 is negatively correlated when using raw ES change
% 2. Section 3 is postiively correlated when using normalized ES change
% both hold when looking at only significant points. section 2 becomes
% postively correlated as well, but can ignore.
% If STDP, would expect section 1 to be postively correlated and section 3
% to be negatively correlated...

ind = {};
change = ESchange; %set this to ESchange or normESchange
keep = Cond;
ind{1} = normISI < 0 & keep & change < 0;%normISI < 0 & ~Bad;
ind{2} = normISI > 0 & normISI < 3 & change > 0 & keep;%normISI < 1 & normISI > 0 & ESchange > 0 & ~Bad;
ind{3} = normISI > 0 & normISI < 4 & change < 0 & keep;%normISI > 0 & normISI < 10 & ESchange < 0 & ~Bad;
ind{4} = ~ind{1} & ~ind{2} & ~ind{3} & keep;

colors = {'r','g','b','m'};

% Aggragate box plots of inhibition change
figure; pvals = [];
for i = 1:4
    hold on; boxplot(IHchange(ind{i}),'colors',colors{i},'position',i,'notch','on','widths',0.5);
    [pvals(i),~] = signrank(IHchange(ind{i}));
end
xlim([0.5,4.5]); xticks(1:4); xticklabels(1:4);
xl = xlim; hold on; plot(xl,[0,0],'k--');

% Scatter plot of evoked spike changes but colored by section
figure; subplot(2,2,1);
for i = 1:4
    scatter(normISI(ind{i} & monkeyJ),change(ind{i} & monkeyJ),45,[colors{i},'^'],'linewidth',1.5);
    hold on; scatter(normISI(ind{i} & monkeyL),change(ind{i} & monkeyL),40,[colors{i},'o'],'linewidth',1.5);
end
xl = xlim; yl = ylim;
hold on; plot([0,0],yl,'k--'); plot(xl,[0,0],'k--');
xlim(xl); ylim(yl);  box off;
xlabel(''); ylabel('Evoked spike % Change'); xlabel('Normalized ISI');

% Scatter plot of inhibition changes but colored by section
subplot(2,2,2);
for i = 1:4
    scatter(normISI(ind{i} & monkeyJ),IHchange(ind{i} & monkeyJ),45,[colors{i},'^'],'linewidth',1.5);
    hold on; scatter(normISI(ind{i} & monkeyL),IHchange(ind{i} & monkeyL),40,[colors{i},'o'],'linewidth',1.5);
end
xl = xlim; yl = ylim;
hold on; plot([0,0],yl,'k--'); plot(xl,[0,0],'k--');
xlim(xl); ylim(yl);  box off;
xlabel(''); ylabel('Inhibition % Change'); xlabel('Normalized ISI');

% Scatter plot of evoked spike changes vs inhibition changes
subplot(2,2,3);
for i = 1:4
    scatter(ESchange(ind{i} & monkeyJ),IHchange(ind{i} & monkeyJ),45,[colors{i},'^'],'linewidth',1.5);
    hold on; scatter(ESchange(ind{i} & monkeyL),IHchange(ind{i} & monkeyL),40,[colors{i},'o'],'linewidth',1.5);
end
xl = xlim; yl = ylim;
hold on; plot([0,0],yl,'k--'); plot(xl,[0,0],'k--');
xlim(xl); ylim(yl);  box off;
xlabel(''); ylabel('Inhibition % Change'); xlabel('Evoked Spike % Change');
title('Raw ES change');
b = []; m = []; p = [];
for i = 1:4
    mdl = fitlm(ESchange(ind{i}),IHchange(ind{i}),'robustopts','huber');
    b(i) = mdl.Coefficients{1,1}; m(i) = mdl.Coefficients{2,1}; 
    p(i) = mdl.Coefficients{2,4};
    hold on; plot(ESchange(Cond), m(i)*ESchange(Cond)+b(i), [colors{i},'--'], 'linewidth', 1.5);
end

% Scatter plot of normalized evoked spike changes vs inhibition changes
subplot(2,2,4);
for i = 1:4
    scatter(normESchange(ind{i} & monkeyJ),IHchange(ind{i} & monkeyJ),45,[colors{i},'^'],'linewidth',1.5);
    hold on; scatter(normESchange(ind{i} & monkeyL),IHchange(ind{i} & monkeyL),40,[colors{i},'o'],'linewidth',1.5);
end
xl = xlim; yl = ylim;
hold on; plot([0,0],yl,'k--'); plot(xl,[0,0],'k--');
xlim(xl); ylim(yl);  box off;
xlabel(''); ylabel('Inhibition % Change'); xlabel('Evoked Spike % Change');
title('Normalized ES change');
b = []; m = []; p = [];
for i = 1:4
    mdl = fitlm(normESchange(ind{i}),IHchange(ind{i}),'robustopts','huber');
    b(i) = mdl.Coefficients{1,1}; m(i) = mdl.Coefficients{2,1}; 
    p(i) = mdl.Coefficients{2,4};
    hold on; plot(normESchange(Cond), m(i)*normESchange(Cond)+b(i), [colors{i},'--'], 'linewidth', 1.5);
end

%% Over time 
%% Evoked spikes over time
% Group 3 decrease diminishes over time

% PRE EPOCH
PreOverTime = {};
bw = 10; 
for i = 1:length(PreES)
    if(isnan(PreES{i})), continue; end
    bins = 0:bw:(PreStim{i}(end));
    Stimcount = histcounts(PreStim{i},bins);
    EScount = histcounts(PreSpikes{i}(PreES{i}),bins);
    PreOverTime{i} = EScount./Stimcount;
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
PreAllZ = (PreAll - nanmean(PreAll)) ./ nanstd(PreAll);

% plot
x = linspace(bw/2,temp*bw-bw/2,temp);
figure; plot(x,PreAllZ,'color',[0.8,0.8,0.8]);
avg = nanmean(PreAllZ,2); stddev = nanstd(PreAllZ,[],2);
hold on; plot(x,avg,'k','linewidth',2);
hold on; plot(x,avg + stddev ,'r');
hold on; plot(x,avg - stddev ,'r');
xlabel('Time (s)'); ylabel('Standard Deviations');
ylim([-4,4]); set(gca,'FontSize',12);

% POST EPOCH
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

PostAllZ = (PostAll - nanmean(PostAll)) ./ nanstd(PostAll);

% Plot by sections
ind = {};
change = ESchange; 
keep = Cond;
ind{1} = normISI < 0 & keep & change < 0;%normISI < 0 & ~Bad;
ind{2} = normISI > 0 & normISI < 3 & change > 0 & keep;%normISI < 1 & normISI > 0 & ESchange > 0 & ~Bad;
ind{3} = normISI > 0 & normISI < 4 & change < 0 & keep;%normISI > 0 & normISI < 10 & ESchange < 0 & ~Bad;
ind{4} = ~ind{1} & ~ind{2} & ~ind{3} & keep;

x = linspace(bw/2,temp*bw-bw/2,temp);
figure;
for i = 1:4
    subplot(1,4,i);
    plot(x,PostAllZ(:,ind{i}),'color',[0.8,0.8,0.8]);
    avg = nanmean(PostAllZ(:,ind{i}),2); stddev = nanstd(PostAllZ(:,ind{i}),[],2);
    hold on; plot(x,avg,'k','linewidth',2);
    hold on; plot(x,avg + stddev ,'r');
    hold on; plot(x,avg - stddev ,'r');
    xlabel('Time (s)'); ylabel('Standard Deviations');
    ylim([-4,4]); set(gca,'FontSize',12);
end

% Box plots and statiscal testing between first and last minute
figure;
for i = 1:4
    subplot(1,4,i);
    first = PostAllZ(1:6,ind{i}); first = first(:);
    last = PostAllZ(end-5:end,ind{i}); last = last(:);
    boxplot([first,last],[zeros(1,length(first)),ones(1,length(last))],'notch','on','symbol','w');
    [p,h] = ranksum(first,last);
    sigstar({[1,2]},p,0,20,0); 
end

%% Inhibition over time
% Group 1 decreases over time

% PRE EPOCH
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
PreAllZ = (PreAll - nanmean(PreAll)) ./ nanstd(PreAll);

% plot
x = linspace(bw/2,temp*bw-bw/2,temp);
figure; plot(x,PreAllZ,'color',[0.8,0.8,0.8]);
avg = nanmean(PreAllZ,2); stddev = nanstd(PreAllZ,[],2);
hold on; plot(x,avg,'k','linewidth',2);
hold on; plot(x,avg + stddev ,'r');
hold on; plot(x,avg - stddev ,'r');
xlabel('Time (s)'); ylabel('Standard Deviations');
ylim([-4,4]); set(gca,'FontSize',12);

% POST EPOCH
PostOverTime = {};
bw = 10;
for i = 1:length(PostInhib)
    bins = 0:bw:(PostStim{i}(end));
    Stimind = discretize(PostStim{i},bins);
    PostOverTime{i} = arrayfun(@(x) nanmedian(PostInhib{i}(Stimind==x)), 1:max(Stimind));
end

% Convert to matrix
temp = min([max(cellfun(@length, PostOverTime)),60]);
PostAll = nan(temp,length(PostOverTime));
for i = 1:length(PostOverTime)
    if(isempty(PostOverTime{i})),continue; end
    inds = 1:min(length(PostOverTime{i}),temp);
    PostAll(inds,i) = PostOverTime{i}(inds);
end

PostAllZ = (PostAll - nanmean(PostAll)) ./ nanstd(PostAll);

% Plot by sections
ind = {};
change = ESchange; 
keep = Cond;
ind{1} = normISI < 0 & keep & change < 0;%normISI < 0 & ~Bad;
ind{2} = normISI > 0 & normISI < 3 & change > 0 & keep;%normISI < 1 & normISI > 0 & ESchange > 0 & ~Bad;
ind{3} = normISI > 0 & normISI < 4 & change < 0 & keep;%normISI > 0 & normISI < 10 & ESchange < 0 & ~Bad;
ind{4} = ~ind{1} & ~ind{2} & ~ind{3} & keep;

x = linspace(bw/2,temp*bw-bw/2,temp);
figure;
for i = 1:4
    subplot(1,4,i);
    plot(x,PostAllZ(:,ind{i}),'color',[0.8,0.8,0.8]);
    avg = nanmean(PostAllZ(:,ind{i}),2); stddev = nanstd(PostAllZ(:,ind{i}),[],2);
    hold on; plot(x,avg,'k','linewidth',2);
    hold on; plot(x,avg + stddev ,'r');
    hold on; plot(x,avg - stddev ,'r');
    xlabel('Time (s)'); ylabel('Standard Deviations');
    ylim([-4,4]); set(gca,'FontSize',12);
end

% Box plots and statiscal testing between first and last minute
figure;
for i = 1:4
    subplot(1,4,i);
    first = PostAllZ(1:6,ind{i}); first = first(:);
    last = PostAllZ(end-5:end,ind{i}); last = last(:);
    boxplot([first,last],[zeros(1,length(first)),ones(1,length(last))],'notch','on','symbol','w');
    [p,h] = ranksum(first,last);
    sigstar({[1,2]},p,0,20,0); 
end

%% First and last minute comparison plots of ES
% Evoked spikes
% POST EPOCH
PreESProb = cellfun(@(x,y) sum(x)/length(y), PreES, PreStim);
PostOverTime = {};
bw = 10; 
for i = 1:length(PostES)
    if(isnan(PostES{i})), continue; end
    bins = 0:bw:(PostStim{i}(end));
    Stimcount = histcounts(PostStim{i},bins);
    EScount = histcounts(PostSpikes{i}(PostES{i}),bins);
    PostOverTime{i} = EScount./Stimcount;
end
% Can't do cellfun cause of empty cells
FirstESProb = [];
LastESProb = [];
for i = 1:length(PostOverTime)
    if(isempty(PostOverTime{i})), continue; end
    FirstESProb(i) = nanmean(PostOverTime{i}(1:6));
    LastESProb(i) = nanmean(PostOverTime{i}(end-5:end));
end
FirstESChange = (FirstESProb-PreESProb)./PreESProb*100;
LastESChange = (LastESProb-PreESProb)./PreESProb*100;

% Plot - Doesn't seem like this is the best way to present the data
% Calculate normalized ISI
normISI = cellfun(@(x,y) x/1000/nanmedian(y), mat2cell(ISI,1,ones(1,length(ISI))), Predt);

subplot(1,2,1); scatter(normISI(~Bad & monkeyJ),FirstESChange(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),FirstESChange(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('First Minute of Post')

legend({'Monkey J','Monkey L'});

subplot(1,2,2); scatter(normISI(~Bad & monkeyJ),LastESChange(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),LastESChange(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Last Minute of Post')

legend({'Monkey J','Monkey L'});

% Plot just the change
ESChangeDiff = LastESChange-FirstESChange;
figure; scatter(normISI(~Bad & monkeyJ),ESChangeDiff(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),ESChangeDiff(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Difference between First and Last minute')

%% First and last minute comparison plots of Inhib
% Inhibition
% PRE EPOCH
PreIHProb = cellfun(@nanmedian, PreInhib);
% POST EPOCH
PostOverTime = {};
bw = 10;
for i = 1:length(PostInhib)
    bins = 0:bw:(PostStim{i}(end));
    Stimind = discretize(PostStim{i},bins);
    PostOverTime{i} = arrayfun(@(x) nanmedian(PostInhib{i}(Stimind==x)), 1:max(Stimind));
end
FirstIHProb = cellfun(@(x) nanmean(x(1:6)),PostOverTime);
LastIHProb = cellfun(@(x) nanmean(x(end-5:end)),PostOverTime);

FirstIHChange = (FirstIHProb-PreIHProb)./PreIHProb*100;
LastIHChange = (LastIHProb-PreIHProb)./PreIHProb*100;

% Plot - Doesn't seem like this is the best way to present the data
% Calculate normalized ISI
normISI = cellfun(@(x,y) x/1000/nanmedian(y), mat2cell(ISI,1,ones(1,length(ISI))), Predt);

subplot(1,2,1); scatter(normISI(~Bad & monkeyJ),FirstIHChange(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),FirstIHChange(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('First Minute of Post')

legend({'Monkey J','Monkey L'});

subplot(1,2,2); scatter(normISI(~Bad & monkeyJ),LastIHChange(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),LastIHChange(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Last Minute of Post')

legend({'Monkey J','Monkey L'});

% Plot just the change
IHChangeDiff = LastIHChange-FirstIHChange;
figure; scatter(normISI(~Bad & monkeyJ),IHChangeDiff(~Bad & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(~Bad & monkeyL),IHChangeDiff(~Bad & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Difference between First and Last minute')







