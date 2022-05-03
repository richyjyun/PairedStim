function PlotOverTime(x,data,comparisontype,ylab,yl,leg,maxcompare)

% Repeated measures ANOVA
tempdata = data';
meas = num2cell(1:size(tempdata,2)); 
meas = cellfun(@num2str,meas,'uniformoutput',false);
labels = cellfun(@(x) ['t',x],meas,'uniformoutput',false);

t = array2table(tempdata,'VariableNames',labels);

rm = fitrm(t,'t1-t60~1');
c = multcompare(rm,'Time','ComparisonType',comparisontype);

% Plot traces and average
a = plot(x,data,'color',[0.8,0.8,0.8]); 
avg = nanmean(data,2); stddev = nanstd(data,[],2);
hold on; 
b = plot(x,avg,'k','linewidth',2);
xlabel('Time (s)'); ylabel(ylab); ylim(yl); 

% Plot significance 
sig = [];
for compare = 1:maxcompare
    inds = find(c{:,1}==compare);
    siginds = find(c{inds,5}< 0.05);
    sig = union(sig,c{inds(siginds),2});
end
for i = 1:length(sig)
    hold on; plot([x(sig(i))-5,x(sig(i))+5],[0.75,0.75],'r','linewidth',2)
end

% Legend
if(leg)
    legend([a(1),b], {'Traces','Average'}, 'box', 'off');
end

set(gca,'FontSize',10); box off;

end