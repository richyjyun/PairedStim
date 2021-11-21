%% Make the session list of all experiments

clear; close all;

user = getenv('username');

metafile = 'C:\Log.xlsx';
opts = detectImportOptions(metafile);
opts.VariableTypes{1} = 'char';
metadata = readtable(metafile,opts);

SLfile = 'C:\SL.mat';

SL = struct([]);

fs = 30000;
for m = 1:size(metadata,1)
    try
        data = metadata(m,:);
        
        if(isempty(data.Date{1})), continue; end
        
        SD = struct();
        SD.Date = data.Date{1};
        
        disp(SD.Date);
        SD.Animal = cell2mat(data.Animal);
        SD.ISI = data.ISI;
        SD.Bad = data.Condition == 0;
        SD.Control = data.Control;
        SD.BadEP = data.BadEP==1;
        
        %% Set path and variables
        path = ['R:\Yun\',SD.Animal,'\Ripple\'];
        datapath = fullfile(path,num2str(SD.Date));
        files = dir(fullfile(datapath,'\PairedPulse*.nev'));
        filename = files.name(1:end-4);
        temp = split(filename,'_');
        temp = split(temp{2},'to');
        
        prechn = str2num(temp{1});
        postchn = str2num(temp{2});
        
        SD.PreChn = prechn;
        SD.PostChn = postchn;
        
        %% Loading stim data
        fprintf('Loading Stim...');
        tic;
        
        stimfile = fullfile(datapath,[filename,'_StimRaw.mat']);
        if(~exist(stimfile))
            Data = LoadRippleData('path',datapath,'filename',filename,'flag',4);
            names = fieldnames(Data);
            for i=1:length(names)
                if(~isempty(Data.(names{i})))
                    eval([names{i} '=Data.' names{i} ]);
                end
            end
            save(stimfile,'stim','stimchannel');
        else
            load(stimfile);
        end
        
        % Set pre and post stim
        preind = find(stimchannel==prechn);
        postind = find(stimchannel==postchn);
        prestim = stim.timestamp(preind,stim.timestamp(preind,:)<min(stim.timestamp(postind,:)));
        prestim = prestim(1:end-1);
        poststim = stim.timestamp(preind,stim.timestamp(preind,:)>max(stim.timestamp(postind,:)));
        condstim = stim.timestamp(postind,:); condstim(isnan(condstim)) = [];
        
        t = toc;
        fprintf('%f seconds\n',t);
        
        %% Load in pre and post test
        fprintf('Evoked Spikes and CCEPs...');
        tic;
        
        datafile = fullfile(datapath,[filename,'_PairedChns.mat']);
        if(exist(datafile))
            load(datafile)
        else
            [prewave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~] = ...
                LoadRippleData_Old('path',datapath,'filename',filename,'flag',1,...
                'times',[prestim(1)-0.5,prestim(end)+0.5],'wavechannel',postchn);
            prestim_shift = prestim-(prestim(1)-0.5);
            
            [postwave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~] = ...
                LoadRippleData_Old('path',datapath,'filename',filename,'flag',1,...
                'times',[poststim(1)-0.5,poststim(end)+0.5],'wavechannel',postchn);
            poststim_shift = poststim-(poststim(1)-0.5);
            
            save(datafile,'prewave','postwave','prestim_shift','poststim_shift','-v7.3');
        end
        
        SD.PreStim = prestim_shift;
        SD.PostStim = poststim_shift;
        
        %% Sort and save spikes
        savefile = fullfile(datapath,[filename,'_PairedSpikes.mat']);
        if exist(savefile)
            load(savefile)
        else
            % Filter
            prefilt = HPF(prewave,fs,1500);
            prefilt = LPF(prefilt,fs,2000);
            postfilt = HPF(postwave,fs,1500);
            postfilt = LPF(postfilt,fs,2000);
            
            spiketimes{1} = SortSpikes(prefilt,fs,SD.PreStim);
            spiketimes{2} = SortSpikes(postfilt,fs,SD.PostStim);
            save(savefile,'prestim_shift','poststim_shift','spiketimes','-v7.3');
        end
        
        SD.PreSpikes = spiketimes{1};
        SD.PostSpikes = spiketimes{2};
        
        %% Save CCEPs
        range = round(-0.2*fs):round(0.3*fs);
        
        trig = round(SD.PreStim*fs);
        [trialinds,bad1] = getTrialinds(trig,range,length(prewave));
        SD.PreEP = mean(prewave(trialinds),2);
        SD.PreEPstd = std(prewave(trialinds),[],2);
        
        trig = round(SD.PostStim*fs);
        [trialinds,bad1] = getTrialinds(trig,range,length(postwave));
        SD.PostEP = mean(postwave(trialinds),2);
        SD.PostEPstd = std(postwave(trialinds),[],2);

        SD.range = range;
        
        % Filtered
        prefilt = bpfilt(prewave,[0.1,500],fs,2);
        postfilt = bpfilt(postwave,[0.1,500],fs,2);
        
        trig = round(SD.PreStim*fs);
        [trialinds,bad1] = getTrialinds(trig,range,length(prefilt));
        PreEP = prefilt(trialinds);
        SD.PreEPFilt = mean(PreEP,2);
        SD.PreEPFiltstd = std(PreEP,[],2);
        
        trig = round(SD.PostStim*fs);
        [trialinds,bad1] = getTrialinds(trig,range,length(postfilt));
        PostEP = postfilt(trialinds);
        SD.PostEPFilt = mean(PostEP,2);
        SD.PostEPFiltstd = std(PostEP,[],2);
        
        SD.SmallerSpike = metadata.SmallerSpike(m);
        
        % Over time of filtered
        bw = 10; bins = 0:bw:max(SD.PreStim);
        inds = discretize(SD.PreStim,bins);
        SD.PreEPTime = arrayfun(@(x) mean(PreEP(:,inds==x),2), 1:max(inds), 'UniformOutput', false);
        SD.PreEPTime = cell2mat(SD.PreEPTime);
        
        bw = 10; bins = 0:bw:max(SD.PostStim);
        inds = discretize(SD.PostStim,bins);
        SD.PostEPTime = arrayfun(@(x) mean(PostEP(:,inds==x),2), 1:max(inds), 'UniformOutput', false);
        SD.PostEPTime = cell2mat(SD.PostEPTime);

        t = toc;
        fprintf('%f seconds\n',t);
        
        %% Load data before and after stim, calculate coherence
        fprintf('Coherence...');
        tic;
        
        allstim = sort(stim.timestamp(:));
        
        if(allstim(1) < 100)
            SD.PreCoh = [];
            SD.PostCoh = [];
            SD.Cohf = [];
        else
            cohfile = fullfile(datapath,[filename,'_Coherence.mat']);
            if(exist(cohfile))
                load(cohfile)
            else
                
                datafile = fullfile(datapath,[filename,'_PairedChnsNoStim.mat']);
                if(exist(datafile))
                    load(datafile)
                else
                    [prewave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~] = ...
                        LoadRippleData_Old('path',datapath,'filename',filename,'flag',1,...
                        'times',[0,floor(min(allstim))],'wavechannel',[prechn,postchn]);
                    
                    [postwave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~] = ...
                        LoadRippleData_Old('path',datapath,'filename',filename,'flag',1,...
                        'times',[ceil(max(allstim))+1,inf],'wavechannel',[prechn,postchn]);
                    
                    save(datafile,'prewave','postwave','-v7.3');
                end
                
                [precoh,f] = mscohere(prewave(:,1),prewave(:,2),10*fs,1*fs,[],fs);
                [postcoh,f] = mscohere(postwave(:,1),postwave(:,2),10*fs,1*fs,[],fs);
                
                save(cohfile,'f','precoh','postcoh','-v7.3');
                
            end
            
            SD.PreCoh = precoh;
            SD.PostCoh = postcoh;
            SD.Cohf = f;
            
        end
        
        t = toc;
        fprintf('%f seconds\n',t);
        
        %% Consolidate
        if(isempty(SL))
            SL = SD;
        else
            SL(end+1) = SD;
        end
    catch
        fprintf('ERROR\n');
    end
end

alldates = extractfield(SL,'Date');
[~,ind] = sort(alldates);
SL = SL(ind);

save(SLfile,'SL','-v7.3');







