function [] = plx2mda(fileName,chansPerTrode,directory)

if nargin == 1
    chansPerTrode = 1;
    directory = pwd;
elseif nargin == 2
    directory = pwd;
end

cd(directory);

if exist(strcat(fileName(1:end-4),'.mat'),'file') ~= 2
    readall(fileName);
    pause(1);
end

load(strcat(fileName(1:end-4),'.mat'),'allts','allwaves',...
    'wfcounts','Freq','nunits1','npw','allad','svStrobed','tsevs',...
    'adfreqs','DateTime','SlowPeakV','adgains','SlowADResBits');

% REORGANIZE SPIKING DATA
temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);
totalUnits = sum(sum(temp));

unitChannel = zeros(totalUnits,1);
temp = cell(totalUnits,1);
temp2 = cell(totalUnits,1);
count = 1;
for ii=1:numChans
   for jj=1:nunits1
       if isempty(allts{jj,Chans(ii)}) == 0
           temp{count} = allts{jj,Chans(ii)};
           unitChannel(count) = ii;
           temp2{count} = allwaves{jj,Chans(ii)};
           count = count+1;
       end
   end
end

allts = temp;
allwaves = temp2;

temp = find(~cellfun(@isempty,allad));
adChans = length(temp);

if adChans == numChans
   check = length(allad{temp(1)});
   
   if check>10
       newad = cell(adChans/chansPerTrode,1);
       newgains = zeros(adChans/chansPerTrode,1);
       newfreqs = zeros(adChans/chansPerTrode,1);
       
       count = 1;
       for ii=1:chansPerTrode:adChans
           newad{count} = allad{temp(ii)};
           newgains(count) = adgains(temp(ii));
           newfreqs(count) = adfreqs(temp(ii));
           count = count+1;
       end
       allad = newad;
       adgains = newgains;
       adfreqs = newfreqs;
   else
       allad = cell(numChans/chansPerTrode,1);
       adgains = adgains(33).*ones(numChans/chansPerTrode,1);
       adfreqs = adfreqs(33).*ones(numChans/chansPerTrode,1);
   end
elseif adChans>0
    check = length(allad{temp(1)});
    if check>10
        newad = cell(numChans/chansPerTrode,1);
        for ii=1:numChans/chansPerTrode
           newad{ii} = allad{temp(1)};
        end
        allad = newad;
        adgains = adgains(33).*ones(numChans/chansPerTrode,1);
        adfreqs = adfreqs(33).*ones(numChans/chansPerTrode,1);
    else
       allad = cell(numChans/chansPerTrode,1);
       adgains = adgains(33).*ones(numChans/chansPerTrode,1);
       adfreqs = adfreqs(33).*ones(numChans/chansPerTrode,1);
    end
else
   allad = cell(numChans/chansPerTrode,1);
   adgains = adgains(33).*ones(numChans/chansPerTrode,1);
   adfreqs = adfreqs(33).*ones(numChans/chansPerTrode,1);
end

numUniqueUnits = totalUnits/chansPerTrode;

newwaves = cell(numUniqueUnits,1);
newts = cell(numUniqueUnits,1);

count = 1;
for ii=1:numUniqueUnits
    times = allts{count};
    waves = allwaves{count};
    waves = waves';
    [numSamples,numEvents] = size(waves);
    fullWaves = zeros(chansPerTrode,numSamples,numEvents);
    
    minicount = 0;
    fullAd = zeros(size(allad{1}),1);
    for jj=1:chansPerTrode
        waves = allwaves{count};
        waves = waves';
        fullWaves(jj,:,:) = waves;
        
        temp = allad{count};
        count = count+1;
    end
    newwaves{ii} = fullWaves;
    newts{ii} = times;
end

allwaves = newwaves;
allts = newts;
clear temp temp2 newwaves totalUnits fullWaves newts;

mkdir(sprintf('%s',fileName(1:end-4)));
cd(sprintf('%s',fileName(1:end-4)));

allEventTimes = cell(numUniqueUnits,1);
for ii=1:numUniqueUnits
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
   
   raw=reshape(cat(2,tmpwaves,zeros(size(tmpwaves))),chansPerTrode,2*npw*size(tmpwaves,3));
   event_times=peak_inds+(0:size(tmpwaves,3)-1)*2*npw;
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

clear ii i j jj q raw event_times tmpwaves maxpeak peak_inds peaks

save(sprintf('%s-mda.mat',fileName(1:end-4)));