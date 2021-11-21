%% Analyzing CCEPs in all paired stimulation experiments. 
% Run cs_AnalyzeStimResponse to run comparisons with evoked spikes and
% inhibition

%% Load data
if(~exist('SL'))
    SLfile = 'C:\SL.mat';
    load(SLfile);
end

%% Extract all variables
ISI = extractfield(SL,'ISI');
Animal = extractfield(SL,'Animal');
BadEP = extractfield(SL,'BadEP');
BadEP = cellfun(@(x) x==1, BadEP);
PreEP = {}; PostEP = {}; Range = {};
PreEPTime = {}; PostEPTime = {};
PreEPstd = {}; PostEPstd = {};
fs = 30000;
filt = 1; % whether we want filted (0.1-500Hz) EPs or raw EPs
for i = 1:length(SL)
    if(filt)
        PreEP{i} = SL(i).PreEPFilt;
        PostEP{i} = SL(i).PostEPFilt;
        PreEPstd{i} = SL(i).PreEPFiltstd;
        PostEPstd{i} = SL(i).PostEPFiltstd;
    else
        PreEP{i} = SL(i).PreEP;
        PostEP{i} = SL(i).PostEP;
        PreEPstd{i} = SL(i).PreEPstd;
        PostEPstd{i} = SL(i).PostEPstd;
    end
    PreEPTime{i} = SL(i).PreEPTime;
    PostEPTime{i} = SL(i).PostEPTime;
    Range{i} = SL(i).range+round(0.0007*fs);
end

% Set indices for analysis
Bad = cell2mat(extractfield(SL,'Bad'));
Control = extractfield(SL,'Control');
Cond = ~(Bad|~isnan(Control));

% Set indices for animals
monkeyJ = cellfun(@(x) strcmp(x,'Jafar'), Animal);
monkeyL = cellfun(@(x) strcmp(x,'Lorde'), Animal);

%% Example figure (81, 85 are good)
x = Range{85}/fs*1000;
blank = find(x>=0,1):find(x>=1.5,1);
y = PreEP{85}; y(blank) = nan;
figure; plot(x,y,'k','linewidth',2);
xlim([-50,100]); xlabel('Time (ms)'); ylabel('Voltage (\muV)');
yl = ylim; 
hold on; fill([0,0,1.5,1.5],[yl,fliplr(yl)],'r',...
    'edgealpha',0,'facealpha',0.3);
box off;

%% Calculate EP magnitudes
% Peak to peak magnitude and corresponding standard deviations for stats
PreEPp2p = []; PostEPp2p = []; 
PreEPp2pstd = []; PostEPp2pstd = [];
PreEPdt = []; PostEPdt = [];
for i = 1:length(PreEP)
    [trough,peak,t_ind,p_ind] = getEPMag(Range{i}/fs,PreEP{i},fs);
    if(isnan(peak))
        PreEPp2p(i) = nan; PreEPdt(i) = nan; PreEPp2pstd(i) = nan;
    else
        PreEPp2p(i) = peak-trough; PreEPdt(i) = (p_ind-t_ind)/fs;
        PreEPp2pstd(i) = sqrt(PreEPstd{i}(t_ind).^2 + PreEPstd{i}(p_ind).^2); 
    end
    
    [trough,peak,t_ind,p_ind] = getEPMag(Range{i}/fs,PostEP{i},fs);
    if(isnan(peak))
        PostEPp2p(i) = nan; PostEPdt(i) = nan; PostEPp2pstd(i) = nan;
    else
        PostEPp2p(i) = peak-trough; PostEPdt(i) = (p_ind-t_ind)/fs;
        PostEPp2pstd(i) = sqrt(PostEPstd{i}(t_ind).^2 + PostEPstd{i}(p_ind).^2);
    end
end

% remove all that could not be properly calculated
BadEP = BadEP | isnan(PreEPp2p) | isnan(PostEPp2p);

% Percent change in EP magnitude
EPp2pChange = (PostEPp2p-PreEPp2p)./PreEPp2p*100;

% Statistics (2 sample t-test with different sample size / variances)
PreN = cellfun(@length, PreStim); PostN = cellfun(@length, PostStim);
PreVar = PreEPp2pstd.^2; PostVar = PostEPp2pstd.^2;
S = sqrt(PreVar./PreN + PostVar./PostN);
tstat = (PreEPp2p-PostEPp2p)./S;
dfnum = (PreVar./PreN + PostVar./PostN).^2;
dfden = ((PreVar./PreN).^2)./(PreN-1) + ((PostVar./PostN).^2)./(PostN-1);
df = dfnum./dfden;

pvalEP = 2*(1-tcdf(abs(tstat),df));

sig = pvalEP<0.05; keep = ~BadEP & Cond;
figure; scatter(normISI(keep & monkeyJ & ~sig),EPp2pChange(keep & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyL & ~sig),EPp2pChange(keep & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyJ & sig),EPp2pChange(keep & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyL & sig),EPp2pChange(keep & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('CCEP Peak to Peak')

legend({'Monkey J','Monkey L'});

% Check if dt from peak / trough of EP is different before and after
% conditioning
EPdtChange = (PostEPdt - PreEPdt) *1000 ;
figure; scatter(normISI(keep & monkeyJ), EPdtChange(keep & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyL), EPdtChange(keep & monkeyL),40,'ko','linewidth',1.5);
ylim([-1,1]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
xlabel('Normalized ISI'); ylabel('% Change');
set(gca,'FontSize',12); title('Change in CCEP dt')

legend({'Monkey J','Monkey L'});

% % RMS 
% % Don't need, peak to peak is enough
% leftlim = cellfun(@(x) find(x/fs>=0.0015,1), Range,'uniformoutput',false);
% rightlim = cellfun(@(x) find(x/fs>=0.05,1), Range,'uniformoutput',false);
% 
% PreEPrms = cellfun(@(x,y,z) rms(x(y:z)), PreEP, leftlim, rightlim);
% PostEPrms = cellfun(@(x,y,z) rms(x(y:z)), PostEP, leftlim, rightlim);
% EPrmsChange = (PostEPrms-PreEPrms)./PreEPrms*100;
% 
% figure; scatter(normISI(keep & monkeyJ & ~sig),EPrmsChange(keep & monkeyJ & ~sig),45,'k^','linewidth',1.5);
% hold on; scatter(normISI(keep & monkeyL & ~sig),EPrmsChange(keep & monkeyL & ~sig),40,'ko','linewidth',1.5);
% hold on; scatter(normISI(keep & monkeyJ & sig),EPrmsChange(keep & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
% hold on; scatter(normISI(keep & monkeyL & sig),EPrmsChange(keep & monkeyL & sig),45,'filled','ko','linewidth',1.5);
% ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
% hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k'); 
% xlabel('Normalized ISI'); ylabel('% Change');
% set(gca,'FontSize',12); title('Peak to Peak')
% 
% legend({'Monkey J','Monkey L'});

%% Compare EP change to ES and IH
% Run cs_AnalyzeEvokedSpikes before this
sig = pvalES<0.05;
figure; subplot(1,2,1); 
scatter(EPp2pChange(keep & monkeyJ & ~sig),ESchange(keep & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyL & ~sig),ESchange(keep & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyJ & sig),ESchange(keep & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyL & sig),ESchange(keep & monkeyL & sig),45,'filled','ko','linewidth',1.5);
xlabel('% Change in CCEPs'); ylabel('% Change in ES');
xlim([-100,100]); ylim([-100,100]); xl = xlim; yl = xlim;
hold on; plot(xl,yl,'k');

% Linear fits - No correlation. The two measures are very different.
mdl = fitlm(EPp2pChange(keep),ESchange(keep));
mdl = fitlm(EPp2pChange(keep & sig),ESchange(keep & sig));

sig = pvalIH<0.05;
subplot(1,2,2); 
scatter(EPp2pChange(keep & monkeyJ & ~sig),IHchange(keep & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyL & ~sig),IHchange(keep & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyJ & sig),IHchange(keep & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(EPp2pChange(keep & monkeyL & sig),IHchange(keep & monkeyL & sig),45,'filled','ko','linewidth',1.5);
xlabel('% Change in CCEPs'); ylabel('% Change in Inhib');
xlim([-100,100]); ylim([-100,100]); xl = xlim; yl = xlim;
hold on; plot(xl,yl,'k');

% Linear fits - No correlation. The two measures are very different.
mdl = fitlm(EPp2pChange(keep),IHchange(keep));
mdl = fitlm(EPp2pChange(keep & sig),IHchange(keep & sig));

legend({'Monkey J','Monkey L'});

%% Over time
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

PreMagZ = (PreMag - nanmean(PreMag,2))./nanstd(PreMag,[],2);

figure; plot(PreMagZ');
hold on; plot(nanmean(PreMagZ),'k','linewidth',2);

% Consolidate Post into matrix
N = max(cellfun(@length, PostEPTimep2p));
PostMag = nan(length(PostEPTimep2p),N);
for i = 1:length(PostEPTimep2p)
    PostMag(i,1:length(PostEPTimep2p{i})) = PostEPTimep2p{i};
end

PostMagZ = (PostMag - nanmean(PostMag,2))./nanstd(PostMag,[],2);

figure;
for i = 1:4
    subplot(1,4,i);
    plot(PostMagZ(ind{i},:)');
    hold on; plot(nanmean(PostMagZ(ind{i},:)),'k','linewidth',2);
end

%% First vs last minute EP mag
First = nanmean(PostMag(:,1:6),2);
Last = cellfun(@(x) nanmean(x(:,end-5:end)), PostEPTimep2p);

FirstChange = (First - PreEPp2p)./PreEPp2p * 100;
LastChange = (Last - PreEPp2p)./PreEPp2p * 100;

keep = ~BadEP & Cond;
figure; subplot(1,2,1); scatter(normISI(keep & monkeyJ),FirstChange(keep & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyL),FirstChange(keep & monkeyL),40,'ko','linewidth',1.5);
xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k');
title('First minute EP mag'); ylabel('% Change'); xlabel('Normalized ISI');

subplot(1,2,2); subplot(1,2,2); scatter(normISI(keep & monkeyJ),LastChange(keep & monkeyJ),45,'k^','linewidth',1.5);
hold on; scatter(normISI(keep & monkeyL),LastChange(keep & monkeyL),40,'ko','linewidth',1.5);
ylim([-100,100]); xlim([-25,25]); xl = xlim; yl = ylim; 
hold on; plot([0,0],yl,'k'); plot([1,1],yl,'k--'); plot(xl,[0,0],'k');
title('Last Minute EP mag'); ylabel('% Change'); xlabel('Normalized ISI');










