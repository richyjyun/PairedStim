%% Analysis of control experiments
% Run cs_AnalyzeStimResponse and cs_AnalyzeCCEPs beforehand

%% Load data
if(~exist('SL'))
    SLfile = 'C:\SL.mat';
    load(SLfile);
end

% labels = {'Long Delay', 'Random Stim', 'Small Current', 'No Stim', 'Pre Stim Only'};

%% Run cs_AnalyzeEvokedSpikes and cs_AnalyzeCCEPs to get data
%% Evoked spikes and inhibition
ControlInd = ~isnan(Control);
sig = pvalES < 0.05;
figure; subplot(1,2,1); scatter(Control(ControlInd & monkeyJ & ~sig),ESchange(ControlInd & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyL & ~sig),ESchange(ControlInd & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyJ & sig),ESchange(ControlInd & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyL & sig),ESchange(ControlInd & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xlim([0.5,max(Control)+0.5]); xl = xlim; yl = ylim; 
hold on; plot(xl,[0,0],'k'); xticks(1:max(Control));
xlabel('Conditions'); ylabel('% Change');
set(gca,'FontSize',12); title('Evoked Spikes')

sig = pvalIH < 0.05;
subplot(1,2,2); scatter(Control(ControlInd & monkeyJ & ~sig),IHchange(ControlInd & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyL & ~sig),IHchange(ControlInd & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyJ & sig),IHchange(ControlInd & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(Control(ControlInd & monkeyL & sig),IHchange(ControlInd & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xlim([0.5,max(Control)+0.5]); xl = xlim; yl = ylim; 
hold on; plot(xl,[0,0],'k'); xticks(1:max(Control));
xlabel('Conditions'); ylabel('% Change');
set(gca,'FontSize',12); title('Inhibition')

legend({'Monkey J','Monkey L'});

%% Evoked spikes vs inhibition
figure;
colors = {'r','g','b','m','c'};
for i = 1:max(Control)
    idx = Control == i;
    hold on; scatter(ESchange(idx & monkeyJ),IHchange(idx & monkeyJ),45,[colors{i},'^'],'linewidth',1.5);
    hold on; scatter(ESchange(idx & monkeyL),IHchange(idx & monkeyL),40,[colors{i},'o'],'linewidth',1.5);
%     mdl(i) = fitlm(ESchange(idx),IHchange(idx));
end
ylabel('Inhib % Change'); xlabel('ES % Change');

%% CCEPs
sig = pvalEP<0.05; keep = ~BadEP & ControlInd;
figure; scatter(Control(keep & monkeyJ & ~sig),EPp2pChange(keep & monkeyJ & ~sig),45,'k^','linewidth',1.5);
hold on; scatter(Control(keep & monkeyL & ~sig),EPp2pChange(keep & monkeyL & ~sig),40,'ko','linewidth',1.5);
hold on; scatter(Control(keep & monkeyJ & sig),EPp2pChange(keep & monkeyJ & sig),50,'filled','k^','linewidth',1.5);
hold on; scatter(Control(keep & monkeyL & sig),EPp2pChange(keep & monkeyL & sig),45,'filled','ko','linewidth',1.5);
ylim([-100,100]); xlim([0.5,max(Control)+0.5]); xl = xlim; yl = ylim; 
hold on; plot(xl,[0,0],'k'); xticks(1:max(Control));
xlabel('Conditions'); ylabel('% Change');
set(gca,'FontSize',12); title('CCEP Peak to Peak')

legend({'Monkey J','Monkey L'});

%% CCEPs vs evoked spikes and inhibition
sig = pvalES<0.05; keep = ~BadEP & ControlInd;
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







