function [] = firingsmda2matoe(directory)
% convert the firings.mda file output from MountainSort back to MATLAB
%  format for data analysis (from Open Ephys format)

% call this code from within the directory that contains all of the
%  firings1.mda, firings2.mda etc. files

%  criterion for inclusion of single units
%     1) <=1% of spikes within a 2ms refractory period
%     2) <0.8 correlation with all other recorded neurons
%     3) average firing rate for the whole session > 0.1 Hz and < 100 Hz

% fileID = fopen('fileName.txt','r');
% fileName1 = fscanf(fileID,'%s');

if nargin<1
    directory = pwd;
end

files = dir('firing*.mda');
numFiles = length(files);

spikeTimes = cell(numFiles,1);
totalUnits = 0;
for ii=1:numFiles
    firingFile = files(ii).name;
    A = readmda(firingFile);
    
    units = unique(A(3,:));
    numUnits = length(units);
    
    spikeTimes{ii} = cell(numUnits,1);
    
    count = 1;
    pseudo_event_times = A(2,:);
    unitcode = A(3,:);
    unitIDs = unique(unitcode);
    nunits = length(unitIDs);
    for kk=1:nunits
        spikeTimes{ii}{count} = pseudo_event_times(unitcode==unitIDs(kk))+1;
        totalUnits = totalUnits+1;
        count = count+1;
    end
end

clearvars -except spikeTimes numFiles directory totalUnits

cd ..

if totalUnits>0
    % convert from pseudo event times to experimental time
    newts = cell(totalUnits,1);
    newwaves = cell(totalUnits,1);
    
    fileName1 = strcat(directory,'-mda.mat');
    load(fileName1)
    
    count = 0;
    for ii=1:numFiles
        spikeTimeTempArray = spikeTimes{ii};
        numUnits = size(spikeTimeTempArray,1);
        
        trueEventTimes = allts{ii};
        trueWaves = allwaves{ii};
        correspondingEventIndices = double(allEventTimes{ii}');
        for jj=1:numUnits
            count = count+1;
            
            pseudo_event_times = unique(spikeTimeTempArray{jj});
            indices = zeros(length(pseudo_event_times),1);
            for kk=1:length(pseudo_event_times)
                [difference,ind] = min(abs(correspondingEventIndices-pseudo_event_times(kk)));
                if abs(difference)<20
                    indices(kk) = ind;
                end
            end
            indices = indices(indices>0);
            newwaves{count} = trueWaves(:,:,indices);
            newts{count} = trueEventTimes(indices);
        end
    end
    clear difference ind trueWaves count;
    
    % include and exclude units
    %  criterion for inclusion:
    %     1) <=1% of spikes within a 2ms refractory period
    %     2) <0.8 correlation with all other recorded neurons
    %     3) average firing rate for the whole session > 0.1 Hz
    
    timeMultiplier = 1000;
    
    pointProcessSpikes = zeros(round(totalTime*timeMultiplier),totalUnits);
    for ii=1:totalUnits
        spikeTimes = max(1,round(newts{ii}.*timeMultiplier));
        for jj=1:length(spikeTimes)
            pointProcessSpikes(spikeTimes(jj),ii) = 1;
        end
    end
    
    refractory_cutoff = 2/1000; % 2ms
    refractory_inclusion = 0.01; % 1%
    spikeHz_cutofflow = 0.1;spikeNum_cutofflow = spikeHz_cutofflow*totalTime; % 0.1 Hz
    spikeHz_cutoffhigh = 100;spikeNum_cutoffhigh = spikeHz_cutoffhigh*totalTime; % 100 Hz
    correlation_inclusion = 0.8; % 0.8 correlation between two neurons throughout recording
    toInclude = ones(totalUnits,1);
    
    for ii=1:totalUnits
        spikeTimes = newts{ii};
        for jj=ii+1:totalUnits
            [r,~] = corrcoef(pointProcessSpikes(:,ii),pointProcessSpikes(:,jj));
            if r(1,2) >= correlation_inclusion
                toInclude(jj) = 0;
                pointProcessSpikes(:,ii) = (pointProcessSpikes(:,ii)+pointProcessSpikes(:,jj))>0;
                pointProcessSpikes(:,jj) = 0;
                temp = newts{jj};
                newts{ii} = unique(round([spikeTimes;temp].*timeMultiplier))./timeMultiplier;
                newts{jj} = 0;
            end
        end
    end
    
    for ii=1:totalUnits
        spikeTimes = newts{ii};
        if length(spikeTimes) < spikeNum_cutofflow || length(spikeTimes) > spikeNum_cutoffhigh
            toInclude(ii) = 0;
        end
        isi = diff([0;spikeTimes]);
        %     figure();subplot(2,1,1);plot(spikeTimes);
        %     subplot(2,1,2);histogram(isi);
        criterion1 = sum(isi<=refractory_cutoff)/length(isi);
%         fprintf('\nProportion refractory violations: %3.2e\n',criterion1);
        if criterion1 > refractory_inclusion || isnan(criterion1)
            toInclude(ii) = 0;
        end
    end
    
    totalUnits = sum(toInclude);
    
    if totalUnits>0
        allts = cell(totalUnits,1);
        allwaves = cell(totalUnits,1);
        meanwaves = cell(totalUnits,3);
        
        inds = find(toInclude==1);

        for ii=1:totalUnits
            allts{ii} = newts{inds(ii)};
            allwaves{ii} = newwaves{inds(ii)};
            
            tempwaves = allwaves{ii};
            meanwaves{ii,2} = squeeze(mean(tempwaves,3));
            
            % run the bootstrap to get 95% confidence intervals on waveform
            alpha = 0.05;
            numIter = 1000;
            [~,numSamples,numEvents] = size(tempwaves);
            result1 = zeros(chansPerTrode,numSamples);
            result2 = zeros(chansPerTrode,numSamples);
            for kk=1:chansPerTrode
                data = squeeze(tempwaves(kk,:,:));
                
                bootstrap = zeros(numIter,numSamples);
                for jj=1:numIter
                    myinds = random('Discrete Uniform',numEvents,[numEvents,1]);
                    temp = data(:,myinds);
                    bootstrap(jj,:) = mean(temp,2);
                end
                Q = quantile(bootstrap,[alpha/2,1-alpha/2],1);
                result1(kk,:) = Q(1,:);
                result2(kk,:) = Q(2,:);
            end
            meanwaves{ii,1} = result1;
            meanwaves{ii,3} = result2;
        end
        
        clear pointProcessSpikes temp temp2 ii jj timeMultiplier newts spikeTimes ...
            trueEventTimes spikeTimeTempArray correspondingEventIndices toInclude r ...
            criterion1 allEventTimes inds isi kk index indices nonemptyad pseudo_event_times ...
            newwaves;
        
        fprintf('\nTotal Units: %d\n',totalUnits);
        
        oldFileName = sprintf('%s.mat',fileName1(1:end-8));
        load(oldFileName,'auxData','eventInfo','events','eventTimes','Fs','lpFs','numChans','lowpassTimes');
        
        newFileName = sprintf('%s-mounsort.mat',fileName1(1:end-8));
        save(newFileName,'allts','allwaves','meanwaves','auxData','eventInfo',...
            'events','eventTimes','Fs','lpFs','numChans','lowpassTimes',...
            'chansPerTrode','totalUnits','totalTime');
    end
    
end
end