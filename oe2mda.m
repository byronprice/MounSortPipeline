function [] = oe2mda(fileName,chansPerTrode,directory)

if nargin == 1
    chansPerTrode = 1;
    directory = pwd;
elseif nargin == 2
    directory = pwd;
end

cd(directory);

load(fileName,'allSpikeData','Fs','lowpassTimes');

totalTime = max(lowpassTimes)-min(lowpassTimes)+10/1000;

% REORGANIZE SPIKING DATA
numChans = size(allSpikeData,1);

allts = cell(numChans,1);
allwaves = cell(numChans,1);

for ii=1:numChans
    allwaves{ii} = allSpikeData{ii,1};
    allts{ii} = allSpikeData{ii,2};
end

numUniqueChans = numChans/chansPerTrode;

newwaves = cell(numUniqueChans,1);
newts = cell(numUniqueChans,1);

count = 1;
for ii=1:numUniqueChans
    times = allts{count};
    waves = allwaves{count};
    waves = waves';
    [numSamples,numEvents] = size(waves);
    fullWaves = zeros(chansPerTrode,numSamples,numEvents);
    
    for jj=1:chansPerTrode
        waves = allwaves{count};
        waves = waves';
        fullWaves(jj,:,:) = waves;
        
        count = count+1;
    end
    newwaves{ii} = fullWaves;
    newts{ii} = times;
end

allwaves = newwaves;
allts = newts;
clear newwaves fullWaves newts times waves;

mkdir(sprintf('%s',fileName(1:end-4)));
cd(sprintf('%s',fileName(1:end-4)));

allEventTimes = cell(numUniqueChans,1);
for ii=1:numUniqueChans
%    waves = allwaves{ii}; 
%    waves = waves'; % number of samples in the snippet by number of events,
%             % if a tetrode, should be 4 by number of samples by number of
%             % events
%    tmpwaves = reshape(waves,[chansPerTrode,size(waves,1),size(waves,2)]);
%    clear waves;
   tmpwaves = allwaves{ii};
   
   [peaks,i] = min(tmpwaves,[],2); % i is the column index of the peak of each channel
   [maxpeak,j] = min(peaks,[],1); % j is the row index of the peak across channel peaks
   
   % use j to index i, find the index of the peak of each waveform across channels
   j = reshape(j,1,size(j,3));
   i = reshape(i,chansPerTrode,size(i,3));
   peak_inds =  zeros(1,size(i,2));
   for q = 1:size(i,2)
       peak_inds(q) = i(j(q),q);
   end
   
   raw=reshape(cat(2,tmpwaves,zeros(size(tmpwaves))),chansPerTrode,2*numSamples*size(tmpwaves,3));
   event_times=peak_inds+(0:size(tmpwaves,3)-1)*2*numSamples;
   event_times = int32(event_times);
   
   allEventTimes{ii} = event_times;
   
%    if ~exist(sprintf('%s.mda',fileName(1:end-4)),'dir')
%        mkdir(sprintf('%s.mda',fileName(1:end-4)))
%    end
%    cd(sprintf('%s',fileName(1:end-4)));
   
   %write mda files
   writemda(event_times,sprintf('event_times.nt%02d.mda',ii),'int32');
   writemda(raw,sprintf('raw.nt%02d.mda',ii),'float64');
   
end

cd(directory);

% clear ii i j jj q raw event_times tmpwaves maxpeak peak_inds peaks

save(sprintf('%s-mda.mat',fileName(1:end-4)),'allts','allwaves','Fs',...
    'allEventTimes','totalTime','chansPerTrode');

end