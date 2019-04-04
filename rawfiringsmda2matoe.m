function [] = rawfiringsmda2matoe(directory)
% convert the rawfire.mda file output from MountainSort back to MATLAB
%  format for data analysis (from Open Ephys format)

% call this code from within the directory that contains the file

%  criterion for inclusion of single units
%     1) <=1% of spikes within a 2ms refractory period
%     2) <0.8 correlation with all other recorded neurons
%     3) average firing rate for the whole session > 0.1 Hz and < 100 Hz

% fileID = fopen('fileName.txt','r');
% fileName1 = fscanf(fileID,'%s');

if nargin<1
    directory = pwd;
end

A = readmda('rawfire.mda');
unitIDs = unique(A(3,:));
totalUnits = length(unitIDs);
spikeTimes = cell(totalUnits,1);

pseudo_event_times = A(2,:);
unitcode = A(3,:);
for kk=1:totalUnits
    spikeTimes{kk} = pseudo_event_times(unitcode==unitIDs(kk))+1;
end

clearvars -except spikeTimes directory totalUnits

cd ..

if totalUnits>0
    fileName1 = strcat(directory,'.mat');
    load(fileName1,'totalTime','startTime','Fs',...
        'events','eventInfo','eventTimes','chansPerTrode','auxData',...
        'lowpassTimes')
    timestamps = (startTime:1/Fs:(startTime+totalTime))';
    % convert from pseudo event times to experimental time
    newts = cell(totalUnits,1);

    for ii=1:totalUnits
        spikeTimeTempArray = spikeTimes{ii};
          
        pseudo_event_times = unique(spikeTimeTempArray);
        indices = zeros(length(pseudo_event_times),1);
        for kk=1:length(pseudo_event_times)
            indices(kk) = timestamps(pseudo_event_times(kk));
        end
        indices = indices(indices>0);
        newts{ii} = indices;
    end
    clear indices pseudo_event_times spikeTimeTempArray;
    
    % include and exclude units
    %  criterion for inclusion:
    %     1) <=1% of spikes within a 2ms refractory period
    %     2) <0.8 correlation with all other recorded neurons
    %     3) average firing rate for the whole session > 0.1 Hz
    
    timeMultiplier = 1000;
    
    pointProcessSpikes = zeros(round(totalTime*timeMultiplier),totalUnits);
    for ii=1:totalUnits
        spikeTimes = max(1,round((newts{ii}-startTime).*timeMultiplier));
        for jj=1:length(spikeTimes)
            pointProcessSpikes(spikeTimes(jj),ii) = pointProcessSpikes(spikeTimes(jj),ii)+1;
        end
    end
    
    refractory_cutoff = 1.5/1000; % 2ms
    refractory_inclusion = 0.02; % 2%
    spikeHz_cutofflow = 0.5;spikeNum_cutofflow = spikeHz_cutofflow*totalTime; % 0.1 Hz
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
                newts{ii} = unique(round(sort([spikeTimes;temp]).*timeMultiplier))./timeMultiplier;
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
        
        inds = find(toInclude==1);

        for ii=1:totalUnits
            allts{ii} = newts{inds(ii)};
        end
        
        fprintf('\nTotal Units: %d\n',totalUnits);
        
        oldFileName = sprintf('%s.mat',fileName1(1:end-4));
        load(oldFileName,'auxData','eventInfo','events','eventTimes');
        
        newFileName = sprintf('%s-mounsortfull.mat',fileName1(1:end-4));
        save(newFileName,'allts','auxData','eventInfo',...
            'events','eventTimes','Fs',...
            'chansPerTrode','totalUnits','totalTime','startTime','lowpassTimes');
    else
        fprintf('\nTotal Units: 0\n'); 
    end
    
else
    fprintf('\nTotal Units: 0\n'); 
end
end